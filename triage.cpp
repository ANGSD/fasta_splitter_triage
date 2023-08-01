//g++ triage.cpp  /home/fvr124/metaDMG-cpp/htslib/libhts.a  -lz -I/home/fvr124/metaDMG-cpp/htslib/ -lpthread -lbz2 -lcurl -llzma

#include <zlib.h>
#include <vector>
#include <map>
#include <cstring>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <sys/stat.h>     // for stat, time_t
#include <htslib/kstring.h>
#include <htslib/bgzf.h>
#include <htslib/faidx.h>


struct cmp_str
{
   bool operator()(char const *a, char const *b) const
   {
      return std::strcmp(a, b) < 0;
   }
};

typedef std::map<int,char *> int2char;
typedef std::map<int,int> int2int;
typedef std::map<int,size_t> int2size_t;
typedef std::map<int,std::vector<int> > int2intvec;
typedef std::map<char*,int,cmp_str> char2int;

void strip(char *line) {
  int at = 0;
  for (unsigned i = 0; i < strlen(line); i++)
    if (line[i] == '\t' || line[i] == '\n')
      continue;
    else
      line[at++] = line[i];
  line[at] = '\0';
}

char *strpop(char **str, char split) {
    char *tok = *str;
    while (**str) {
        if (**str != split)
            (*str)++;
        else {
            **str = '\0';
            (*str)++;
            break;
        }
    }
    return tok;
}

void parse_nodes(const char *fname, int2char &rank, int2int &parent, int2intvec &child, int dochild) {
  fprintf(stderr,"\t-> Parsing: \'%s\'\n",fname);
    gzFile gz = Z_NULL;
    gz = gzopen(fname, "rb");
    if (gz == Z_NULL) {
        fprintf(stderr, "\t-> Problems opening file: \'%s\'\n", fname);
        exit(0);
    }
    char buf[4096];
    char **toks = new char *[5];

    while (gzgets(gz, buf, 4096)) {
        strip(buf);  // fprintf(stderr,"buf:%s\n",buf);
        char *saveptr = buf;
        toks[0] = strpop(&saveptr, '|');
        toks[1] = strpop(&saveptr, '|');
        toks[2] = strpop(&saveptr, '|');
        for (int i = 0; 0 && i < 3; i++)
            fprintf(stderr, "%d):\'%s\'\n", i, toks[i]);

        int2int::iterator it = parent.find(atoi(toks[0]));
        if (it != parent.end())
            fprintf(stderr, "\t->[%s] duplicate name(column0): %s\n", fname, toks[0]);
        else {
            int key = atoi(toks[0]);
            int val = atoi(toks[1]);
            parent[key] = val;
            rank[key] = strdup(toks[2]);
            if (dochild) {
                if (key == val)  // catch 1 <-> 1
                    continue;
                int2intvec::iterator it2 = child.find(val);
                if (it2 == child.end()) {
                    std::vector<int> tmp;
                    tmp.push_back(key);
                    child[val] = tmp;
                } else
                    it2->second.push_back(key);
            }
        }
    }
    fprintf(stderr, "\t-> Number of unique names (column1): %lu from file: %s parent.size():%lu child.size():%lu\n", rank.size(), fname, parent.size(), child.size());
    gzclose(gz);
    delete[] toks;
    fprintf(stderr,"\t-> Done reading file: \'%s\'\n",fname);
}

//this function parses the taxid and the genome size, if the genomesize is NA it will be set to some random value between 1e3 and 1e5
//this function assumed the input file follows antonios custom format
int2char taxid2all;

int2size_t parse_meta(const char *fname){
  fprintf(stderr,"\t-> Parsing file: \'%s\'\n",fname);
  int2size_t retval;
  gzFile gz = Z_NULL;
  
  gz = gzopen(fname, "rb");
  if (gz == Z_NULL) {
    fprintf(stderr, "\t-> Problems opening file: \'%s\'\n", fname);
    exit(0);
  }
  char buf[4096];
  char **toks = new char *[10];//tokanize the first 10columns
  gzgets(gz, buf, 4096);  //skip header
  
  while (gzgets(gz, buf, 4096)) {
    char *wholeline = strdup(buf);
    toks[0] = strtok(buf,"\t\n");
    for(int at=1;at<10;at++)
      toks[at] = strtok(NULL,"\t\n");

    int taxid = atoi(toks[6]);
    int2char::iterator it = taxid2all.find(taxid);
    assert(it==taxid2all.end());
    
    taxid2all[taxid] = wholeline;
    if(retval.find(taxid)!=retval.end()){
      fprintf(stderr,"taxid exists: %d %s buf: %s\n",taxid,toks[6],buf);
      for(int i=0;i<10;i++)
	fprintf(stderr,"%d) %s\n",i,toks[i]);
      exit(0);
    }
    size_t size =  atol(toks[8]);
    if(strcmp(toks[8],"NA")==0)
      //      size=0;
    //      size = drand48()*1e6+1000; 

    retval[taxid] = size;
  }
#if 0
  for(auto x: taxid2all) fprintf(stderr,"first.%d last: %s",x.first,x.second);
#endif
  //  exit(0);
  fprintf(stderr,"\t-> Done parsing file: \'%s\', file contains: %lu\n",fname,retval.size());
  return retval;
}


size_t getsum_of_subtree(int2size_t &retmap, int2intvec &child, int taxid) {

    int2size_t::iterator it = retmap.find(taxid);
    if (it != retmap.end())//if we are at leafnode or we have computed the sum of subtree
      return it->second;
    
    size_t ret = 0;
    if (child.size() > 0) {
      int2intvec::iterator it2 = child.find(taxid);
      if (it2 != child.end()) {//loop over subtrees if they exists -1 indicates we should skip that subtree
	std::vector<int> &avec = it2->second;
	for (unsigned i = 0; i < avec.size(); i++) {
	  if(avec[i]!=-1)
	    ret += getsum_of_subtree(retmap, child, avec[i]);
	}
      }
    }
    
    retmap[taxid] = ret;
    return ret;
}

size_t how_many_subnodes(int2intvec &child,int taxid){
  assert(taxid!=-1);

  int2intvec::iterator it = child.find(taxid);
  if(it==child.end())
    return 1;
  std::vector<int> &avec = it->second;
  if(avec.size()==0) //I dont think this can happen
    return 1;
  
  size_t nsum =1;
  for(int i=0;i<avec.size();i++)
    if(avec[i]!=-1){
      nsum += how_many_subnodes(child,avec[i]);
    }
  return nsum;
}


void print_a_subtree(FILE *fp, int2intvec &child,int taxid, int2char &everything){
  //fprintf(stderr,"taxid: %d\n",taxid);
  assert(taxid!=-1);

  int2intvec::iterator it = child.find(taxid);
  if(it==child.end()){
    int2char::iterator it2=everything.find(taxid);
    //    fprintf(stderr,"Im a node with no child: %d\n",taxid);
    if(it2!=everything.end())
      fprintf(fp,"%s",it2->second);
    
    return;
  }
  std::vector<int> &avec = it->second;
  if(avec.size()==0){ //I dont think this can happen
    fprintf(stderr,"I dont think this can happen \n");
    return;
  }
  
  for(int i=0;i<avec.size();i++){
    //fprintf(stderr,"taxid: %d avec[%d]: %d\n",taxid,i,avec[i]);
    if(avec[i]!=-1)
      print_a_subtree(fp,child,avec[i],everything);
  }

}

int print_node(FILE *fp,int2int &up,int2intvec &down,int taxid){
  fprintf(fp,"node:%d ",taxid);
  int2intvec::iterator it2 = down.find(taxid);
  if(it2==down.end())
    return fprintf(fp," No childs\n");
  std::vector<int> &avec = it2->second;
  for(int i=0;i<avec.size();i++)
    fprintf(fp,"(%d)",avec[i]);
  fprintf(fp,"\n");
  for(int i=0;i<avec.size();i++)
    if(avec[i]!=-1)
      print_node(fp,up,down,avec[i]);
  return 0;
}

size_t getval(int2size_t &nucsize,int taxid){
  assert(taxid!=-1);
  int2size_t::iterator it = nucsize.find(taxid);
  assert(it!=nucsize.end());
  return it->second;
}


void getdiffval(int2size_t &nucsize,int taxid,int2intvec &down,std::map<double,int> &mymap,double target){
  //  fprintf(stderr,"[%s] taxid: %d\n",__FUNCTION__,taxid);
  assert(taxid!=-1);assert(taxid!=0);
  
  double mydiff = fabs(getval(nucsize,taxid)-target);

  mymap[mydiff] = taxid;

  int2intvec::iterator it2 = down.find(taxid);
  
  if(it2!=down.end()){
    std::vector<int> &avec = it2->second;
    //fprintf(stderr,"first: %d avec.size(): %lu\n",it2->first,avec.size());
    for(int i=0;i<avec.size();i++)
      if(avec[i]!=-1)
	getdiffval(nucsize,avec[i],down,mymap,target);
  }
}

int fexists(const char *str) {  ///@param str Filename given as a string.
    fprintf(stderr, "\t-> Checking if exits: \'%s\'\n", str);
    struct stat buffer;
    return (stat(str, &buffer) == 0);  /// @return Function returns 1 if file exists.
}

int fexists2(const char *str1, const char *str2) {
    unsigned tmp_l = strlen(str1) + strlen(str2) + 5;
    char tmp[tmp_l];
    snprintf(tmp, tmp_l, "%s%s", str1, str2);
    return fexists(tmp);
}

BGZF *getbgzf(const char *str1, const char *mode, int nthreads) {
    BGZF *fp = NULL;
    fp = bgzf_open(str1, mode);
    fprintf(stderr, "\t-> opening file: \'%s\' mode: \'%s\'\n", str1, mode);
    if (fp == NULL) {
        fprintf(stderr, "\t-> Problem opening file: \"%s\"\n", str1);
        exit(0);
    }
    if (nthreads > 1) {
        fprintf(stderr, "\t-> Setting threads to: %d \n", nthreads);
        bgzf_mt(fp, nthreads, 64);
    }
    return fp;
}

BGZF *getbgzf2(const char *str1, const char *str2, const char *mode, int nthreads) {
    unsigned tmp_l = strlen(str1) + strlen(str2) + 5;
    char tmp[tmp_l];
    snprintf(tmp, tmp_l, "%s%s", str1, str2);
    return getbgzf(tmp, mode, nthreads);
}


int SIG_COND = 1;
char2int acc2taxid(char *acc2taxid_flist ) {
  fprintf(stderr, "\t-> Starting to extract (acc->taxid) from file list: \'%s\'\n", acc2taxid_flist);
  fflush(stderr);
  int dodump = !fexists2(acc2taxid_flist, ".bin");
  
  fprintf(stderr, "\t-> Checking if bimnary file exists. dodump=%d \n", dodump);
  
  time_t t = time(NULL);
  BGZF *fp = NULL;
  kstring_t *kstr = (kstring_t *)malloc(sizeof(kstring_t));
  kstr->l = kstr->m = 0;
  kstr->s = NULL;
  
  if (dodump)
    fp = getbgzf2(acc2taxid_flist, ".bin", "wb", 4);
  else
    fp = getbgzf2(acc2taxid_flist,".bin", "rb", 4);

   char2int am;

   if   (dodump) {
        char buf[4096];
        int at = 0;
        char buf2[4096];
       
        BGZF *fp2 = getbgzf(acc2taxid_flist, "rb", 2);
	kstr->l = 0;
        while (SIG_COND && bgzf_getline(fp2, '\n', kstr)) {
            if (kstr->l == 0)
                break;
            // fprintf(stderr,"at: %d = '%s\'\n",at,kstr->s);
	    if (isatty(fileno(stderr)))
	      fprintf(stderr, "\r\t-> At linenr: %d in \'%s\' entry is: %s \n", at, acc2taxid_flist,kstr->s);
	    BGZF *fp3 = getbgzf(kstr->s, "rb", 2);
	    bgzf_getline(fp3, '\n', kstr);//<- skip header
	    kstr->l = 0;
	    while (SIG_COND && bgzf_getline(fp3, '\n', kstr)) {
	      //	      fprintf(stderr,"streing: %s\n",kstr->s);
	      if (kstr->l == 0)
                break;
	    
	      char *tok = strtok(kstr->s, "\t\n ");
	      char *key = strdup(strtok(NULL, "\t\n "));
	      tok = strtok(NULL, "\t\n ");
	      int val = atoi(tok);
	      //fprintf(stderr,"key: %s val: %d\n",key,val);
	      char2int::iterator it = am.find(key);
	      if (it != am.end())
                fprintf(stderr, "\t-> Duplicate entries found \'%s\' it->first: %s it->second: %d\n", key,it->first,it->second);
	      am[key] = val;
	      // kstr->l = 0;
	    }
	    bgzf_close(fp3);
        }
        bgzf_close(fp2);
	kstr->l=0;
	for(auto it=am.begin();it!=am.end();it++){
	  ksprintf(kstr,"%s\t%d\n",it->first,it->second);
	  if(kstr->l>10000000){
	    bgzf_write(fp,kstr->s,kstr->l);
	    kstr->l =0;
	  }
	}
	bgzf_write(fp,kstr->s,kstr->l);
	
    } else {
        char *key;
	int val;
        while (SIG_COND && bgzf_getline(fp, '\n', kstr)) {
	    if (kstr->l == 0)
	      break;
	    key = strtok(kstr->s,"\t\n ");
	    val = atoi(strtok(NULL,"\t\n "));
	    am[strdup(key)] = val;
        }
    }

    bgzf_close(fp);
    fprintf(stderr, "\t-> Number of entries to use from accesion to taxid: %lu, time taken: %.2f sec\n", am.size(), (float)(time(NULL) - t));
    return am;
}



int get_closest(double target,int2size_t &nucsize, int2intvec &down){
  std::map<double,int> mymap;
  getdiffval(nucsize,1,down,mymap,target);
  std::map<double,int>::iterator it=mymap.begin();
  return it->second;
}


int *splitdb(int nk,int2int &up,int2intvec &down,int2size_t &nucsize){
  size_t total_sum_bp = nucsize.find(1)->second;
  size_t total_sum_nodes = how_many_subnodes(down,1);
  double target = total_sum_bp/nk;
  fprintf(stderr,"\t->[%s] total sum of tree: %lu target for subtrees: %f total_nodes: %lu\n",__FUNCTION__,total_sum_bp,target,total_sum_nodes);
  int *trees = new int[nk];

  for(int k=0;k<nk-1;k++){//loop over subtrees
    fprintf(stderr,"\t->------------ k:%d nk:%d howmany_nodes: %d totalbp: %lu\n",k,nk,how_many_subnodes(down,1),getval(nucsize,1));
    //print_node(stderr,up,down,1);
  
    trees[k] = get_closest(target,nucsize,down);
    size_t subtreebp = nucsize.find(trees[k])->second;
    fprintf(stderr,"\t-> Node to pick: %d how_many_subnods: %lu how_many_in_bp: %lu \n",trees[k],how_many_subnodes(down,trees[k]),subtreebp);


 
    int2int::iterator it = up.find(trees[k]);
    int upnode = it->second;
    fprintf(stderr,"\t-> removing edge: (%d<->%d)\n",it->first,upnode);
    it->second = -1;

    //we need to remove the edge, it is multificating, so we should loop over all child to ensure we remove the proper
    int2intvec::iterator it2 = down.find(upnode);assert(down.end()!=it2);
    std::vector<int> &avec = it2->second;
#if 0
    for(int i=0 ; i<avec.size() ; i++)
      fprintf(stderr,"avec[%d]: %d\n",i,avec[i]);
#endif
    for(int v=0 ; v<avec.size() ; v++){
      if(avec[v]==trees[k]){
	avec[v] = -1;
      }
    }
 
        
    //now loop up through tree and subract all nodes with the value of the subtree we have removed
    it = up.find(upnode);
    //it->first is taxid for the current node where we should subtract the number of bp
    //it->second is parental node
    int lasttaxid = -1;//initialize //stupid way of dealing with breaking from loop at root
    while(1){
      if(lasttaxid==it->first)
	break;

      lasttaxid = it->first;
      fprintf(stderr,"\t-> Looping up tree (%d,%d)\n",it->first,it->second);
      int2size_t::iterator it3 = nucsize.find(it->first);
      assert(it3!=nucsize.end()&&it3->second!=-1&&it3->first!=-1);

      size_t oldsize = it3->second; 
      it3->second -= subtreebp;
      size_t newsize = it3->second;
      assert(it->first==it3->first);
      fprintf(stderr,"\t-> Updating value of taxid(it3->first): %d oldsize: %lu newsize: %lu\n",it->first,oldsize,newsize);
      
      it = up.find(it->second);
    }
  }
  fprintf(stderr,"Done pruning n-1 trees, the remaining tree is the the final subtree which will have the original root as root\n");
  trees[nk-1] = 1;
  size_t subtreebp = nucsize.find(trees[nk-1])->second;
  fprintf(stderr,"\t-> Node to pick: %d how_many_subnods: %lu how_many_in_bp: %lu \n",trees[nk-1],how_many_subnodes(down,trees[nk-1]),subtreebp);


  return trees;
}

//this assumes that that could be different taxids in faifile, like nt
int main_fai_extension1(kstring_t *kstr,gzFile afai,char *afai_fname,char2int &acc2taxid,gzFile fpout){
  char buf[1024];
  char tmp[1024];
  while(gzgets(afai,buf,1024)){
    if(buf[0]=='#')
      continue;
    buf[strlen(buf) - 1] = '\0'; 
    strcpy(tmp,buf);
    char *acc = strtok(buf,"\t \n");
    char2int::iterator it = acc2taxid.find(acc);
    int taxid = -1;
    if(it!=acc2taxid.end())
      taxid = it->second;
    
    ksprintf(kstr,"%s\t%s\t%d\n",afai_fname,tmp,taxid);
    if(1||kstr->l>10000000){
      gzwrite(fpout,kstr->s,kstr->l);
      kstr->l = 0;
    }
  }
  gzwrite(fpout,kstr->s,kstr->l);
  kstr->l = 0;
  return 0;
}

//return taxid second parameter will contain sum of perchrom lengths
int get_total_length(faidx_t *myfai,size_t &total,char2int &acc2taxid){
  int return_taxid = -1;
  total = 0;
  for(int i=0;i<faidx_nseq(myfai);i++){
    const char *namename = faidx_iseq(myfai,i);
    char *name = strdup(namename);//this should be fixed
    char2int::iterator it = acc2taxid.find(name);
    if(it==acc2taxid.end()){
      fprintf(stderr,"\t-> Problem finding name: %s in the acc2taxid map\n",name);
      return -1;
    }
    if(return_taxid==-1)
      return_taxid = it->second;
    else if(return_taxid!=it->second){
      fprintf(stderr,"\t-> Problem with mismatch of taxid for name: %s taxids %d,%d\n",return_taxid,it->second);
      return -1;
    }
    total += faidx_seq_len64(myfai,name);
    free(name);
  }
  return return_taxid;
}


//this assumes that that could be different taxids in faifile, like wgs
int main_fai_extension2(kstring_t *kstr,gzFile afai,char *afai_fname,char2int &acc2taxid,gzFile fpout){
  char buf[1024];
  size_t total = 0;
  while(gzgets(afai,buf,1024)){
    if(buf[0]=='#')
      continue;
    //    fprintf(stderr,"buf: %s\n",buf);
    buf[strlen(buf) - 1-4] = '\0';//last four is for removing .fai so we supply with fasta 
    //    fprintf(stderr,"buf: %s\n",buf);
    faidx_t *myfai = fai_load(buf);
    int taxid = get_total_length(myfai,total,acc2taxid);
    if(taxid==-1){
      fprintf(stderr,"\t-> Error entry: %s\n",buf);
      exit(1);
    }
    ksprintf(kstr,"%s\t%d\t%lu\n",buf,taxid,total);
    //    fprintf(stderr,"%s\t%d\t%lu\n",buf,taxid,total);
    if(1||kstr->l>10000000){
      gzwrite(fpout,kstr->s,kstr->l);
      kstr->l = 0;
    }
    fai_destroy(myfai);
  }
  gzwrite(fpout,kstr->s,kstr->l);
  kstr->l = 0;
  return 0;
}


int main_fai_extension(const char *meta_file,const char *outname,char *acc2taxid_flist,int faitype){
  //  fprintf(stderr,"meta_file: %s outname: %s acc3taxid: %s faitype: %d\n",meta_file,outname,acc2taxid_flist,faitype);
  char2int ass2tax=acc2taxid(acc2taxid_flist);

  gzFile fp = Z_NULL;
  if(((fp=gzopen(meta_file,"rb")))==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s for reading\n",meta_file);
    exit(1);
  }
  fprintf(stderr,"\t-> Opened file: %s for reading\n",meta_file);
  gzFile fpout = Z_NULL;
  if(((fpout=gzopen(outname,"wb")))==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s for writing\n",outname);
    exit(1);
  }
  fprintf(stderr,"\t-> Opened file: %s for writing\n",outname);
  
  char buf[1024];
  kstring_t *kstr = new kstring_t;
  kstr->s = NULL;kstr->l = kstr->m =0;

  while(gzgets(fp,buf,2014)){
    if(buf[0]=='#')
      continue;
    char *tok = strtok(buf,"\n");
    //    fprintf(stderr,"tok: %s\n",tok);
    gzFile fp2 = Z_NULL;
    if(((fp2=gzopen(tok,"rb")))==Z_NULL){
      fprintf(stderr,"\t-> Problem opening file: %s for reading\n",meta_file);
      exit(1);
    }
    if(faitype==1)
      main_fai_extension1(kstr,fp2,tok,ass2tax,fpout);
    else if(faitype ==2)
      main_fai_extension2(kstr,fp2,tok,ass2tax,fpout);
    gzclose(fp2);
  }

  gzclose(fpout);
  gzclose(fp);
  for(auto x: ass2tax) free(x.first);

  return 0;
}



int main(int argc,char **argv){
  const char *node_file = "/projects/lundbeck/scratch/taxDB/v6/taxonomy/nodes.dmp";
  const char *meta_file = "/projects/lundbeck/scratch/taxDB/v6/metadata/taxdb-genome_stats-broad-v6.tsv.gz";
  const char *outname = "outname";
  char *acc2taxid_flist = "acc2taxid_flist";
  int how_many_chunks = 8;
  int makefai = 0;
  for(int at = 1;at<argc;at++){
    if(strcasecmp(argv[at],"-node_file")==0)
      node_file = strdup(argv[at+1]);
    else if(strcasecmp(argv[at],"-meta_file")==0)
      meta_file = strdup(argv[at+1]);
    else if(strcasecmp(argv[at],"-acc2taxid_flist")==0)
      acc2taxid_flist = strdup(argv[at+1]);
    else if(strcasecmp(argv[at],"-nchunks")==0)
      how_many_chunks = atoi(argv[at+1]);
    else if(strcasecmp(argv[at],"-outname")==0)
      outname = strdup(argv[at+1]);
    else if(strcasecmp(argv[at],"makefai")==0){
      makefai = atoi(argv[at+1]);
    }
    at++;;
      
  }
  fprintf(stderr,"\t 1) ./program -node_file filename.txt -meta_file filenames.txt -nchunks integer -acc2taxid_flist file.list\n\t 2) ./program makefai -meta_file filenames.txt -acc2taxid_flist file.list\n\t-> -node_file: \'%s\'\n\t-> -meta_file: \'%s\'\n\t-> -nchunks: %d\n\t-> -outname: %s\n\t-> -acc2taxid_flist: %s\n\t-> makefai: %d\n",node_file,meta_file,how_many_chunks,outname,acc2taxid_flist,makefai);
  if(makefai){
    return main_fai_extension(meta_file,outname,acc2taxid_flist,makefai);
  }
  
  int2char taxid_rank;
  int2int taxid_parent;
  int2intvec taxid_childs;

  parse_nodes(node_file,taxid_rank,taxid_parent,taxid_childs,1);
  //  print_node(stderr,taxid_parent,taxid_childs,1);
  fprintf(stderr,"\t-> howmany nodes: %d taxid_parents.size(): %lu, these values should be identical\n",how_many_subnodes(taxid_childs,1),taxid_parent.size());
  assert(how_many_subnodes(taxid_childs,1)==taxid_parent.size());

#if 0
  for(int2int::iterator it=taxid_parent.begin();it!=taxid_parent.end();it++)
    fprintf(stderr,"\t%d) -> %d\n",it->first,it->second);
#endif
  char2int ass2tax=acc2taxid(acc2taxid_flist);

  return 0;

  
  int2size_t taxid_genome_size =parse_meta(meta_file);
  int2size_t::iterator it = taxid_genome_size.begin();
  fprintf(stderr,"\t-> presize: %lu key: %d val:%lu\n",taxid_genome_size.size(),it->first,it->second);

  getsum_of_subtree(taxid_genome_size, taxid_childs, 1);
  fprintf(stderr,"\t-> postsize: %lu \n",taxid_genome_size.size());
#if 0
  for(int2size_t::iterator it=taxid_genome_size.begin();it!=taxid_genome_size.end();it++)
      if(it->second!=0)
      fprintf(stderr,"\t%d) -> %lu\n",it->first,it->second);
#endif
  //  int *subtrees = splitdb(how_many_chunks,taxid_parent,taxid_childs,taxid_genome_size);
  //subtrees contains the taxids of the "root" of each detached subtree
  for(int i=0;i<how_many_chunks;i++){
    char onam[1024];
    snprintf(onam,1024,"%s.%d",outname,i);
    FILE *fp = fopen(onam,"wb");
    fprintf(stderr,"\t-> Writing file: %s\n",onam);
    print_a_subtree(fp,taxid_childs,1,taxid2all);
    fclose(fp);
    
  }
  
  return 0;
}
