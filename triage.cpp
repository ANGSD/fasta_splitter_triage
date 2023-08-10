//g++ triage.cpp  /home/fvr124/metaDMG-cpp/htslib/libhts.a  -lz -I/home/fvr124/metaDMG-cpp/htslib/ -lpthread -lbz2 -lcurl -llzma
//./a.out -wgs wgs.fix -seqs seqs.fix -node_file nodes_20230719.dmp.gz -nchunks 3 -nrep 11 

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
#include <libgen.h> //for basename




/*



#define MAXDBSIZE;





struct DataBase {
	int DBsize;
	int taxa[MAXDBSIZE]; //you are probably going to set it up differently anyway so I just did this for simplicity. Should really be done with a linked list.
}; //assume staxa are integer valued. You can change it to a set or linked list of character strings.

getlev//finds all the leaf nodes that are descendants of node 'node'
int number_of_leaves_descending_from_node(int node)
{

	int i, k, nn = 0;

	k =number_of_direct_children(node);//this is simply the number of direct children of the node. 
	for (i=0; i<k; i++)
		nn = nn + number_of_leaves_descending_from_node(child node i);
	return nn;
}

//returns a DataBase of size mydbsize for a clade with root node 'node'
//This functions should be called with the root nodes of the different databases in the same order as they were found
struct DataBase Fetch_me_B(int mydbsize, int node)//here I assume nodes are integer valued
	{
	struct Database DBloc, DBret;
	int i, j, v, k, numG;
	
	if (mydbsize==0) DBloc->DBsize = 0;
	else if (node is a leafnode) 
		{
		if (node has no genome) {printf("Encountered leafnode without genome"); exit(-1);}
		else {
			DBloc->DBsize = 1;
			DBloc->taxa[0]=genome in node;
		}

	else {
		k =number_of_direct_children(node);//this is simply the number of direct children of the node. Function depends on how the tree is set up. See 'main()'
		numG = number_of_leaves_descending_from_node(node);
		if (numG < mydbsize) {printf("Cannot find %i genomes\n",mydbsize); exit(-1);}
		DBloc->DBsize = 0;
		for (i=0; i<k; i++)
			{
			v = ceil((double)(mydbsize-DBloc->DBsize)/(double)(k-i)) //how many it should find in the next clade to be on track to find all of them
			numG = numG - number_of_leaves_descending_from_node(child node i); //how many genomes that now maximally can be found after excluding child node i
			while (numG+v < mydbsize-mDBloc->DBsize) v++; //checking that there are enough genomes left and otherwise incrementing
			DBret = Fetch_me_B(v, child node i)  //we could randomize the order of the children?
			for (j=0; j<DBret->DBsize; j++)
				DBloc->taxa[j+DBloc->DBsize] = DBret->taxa[j];
			DBloc->DBsize = DBloc->DBsize + DBret->DBsize;
			}
		}
	return DBloc
	}

//run seperately for euks, prokaryotes and organelle DNA
//We'll do some trial and error to find a good value of c
void main()
{

struct Database Our_database;
int i, c; //number of sequences we want in each database

for (i=0; i<number of databases to be found; i++){ //this should loop over databases in the same order as they were found in the previous program
	Our_database = Fetch_me_B(c, rootnode of database i);
	print the database to file;
	detach rootnode of database i from tree;
	}
}

Some pseudocode.

Objective: to select c phylogenetically diverse genomes for a database

B = {S, k}, k \in integers, S is a set of species (thinking of it as a struct in C, where the set of species would be a linked list or a set in C++. For simplicity of code I assume you can add these together elementwise.)

Fetch_me_B(j, node) {
        if (j==0) return {0, Empty}
        if (node is a leafnode) 
                {
                if (node has no genome) exit with anticipated error. //All leafnodes should have genomes. Otherwise the algorithm needs to be modified slightly
                else return {1, Genome in node}
                }
        else {
                k =number_of_direct_children(node)
                numG = 0
                for i=0,…,k-1
                        numG = numG + number_of_leaves_descending_from_node(child node i)
                myB = {0, Empty} 
                for i=0,…,k-1
                        v = ceiling((j-myB.S)/(k-i)) //how many it should find in the next clade to be on track to find all of them
                        while (numG -v > j-myB.S) v++; //checking that there are enough genomes left and otherwise incrementing
                        myB = myB + Fetch_me_B(v, child number i)  //we could randomize the order of the children
                        }
                }
        return myB
        }
}

 */


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



char2int getlevels() {
  const char *names[49] = {"superkingdom", "domain", "lineage", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "clade", "cohort", "subcohort", "superorder", "order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "infratribe", "genus", "subgenus", "section", "series", "subseries", "subsection", "species", "species group", "species subgroup", "subspecies", "varietas", "morph", "subvariety", "forma", "forma specialis", "biotype", "genotype", "isolate", "pathogroup", "serogroup", "serotype", "strain","no rank"};
  int values[49] = {36, 37, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1};
  
    char2int c2i;
    for (int i = 0; i < 49; i++)
        c2i[strdup(names[i])] = values[i];

    fprintf(stderr, "\t-> Number of entries with level information: %lu \n", c2i.size());
    return c2i;
}

void parse_nodes(const char *fname, int2int &rank, int2int &parent, int2intvec &child, int dochild) {
  fprintf(stderr,"\t-> Parsing: \'%s\'\n",fname);

  char2int c2i = getlevels();
  
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
	    char2int::iterator it = c2i.find(toks[2]);
	    if(it==c2i.end()){
	      fprintf(stderr,"\t-> Problem with rank in file: %s key: %d val: %d rank: %s\n",fname,key,val,toks[2]);
	      rank[key] = 37;
	    }else
	      rank[key] = it->second;
	    
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

size_t getval(int2size_t &nucsize,int taxid){
  //  fprintf(stderr,"taxid: %d\n",taxid);
  assert(taxid!=-1);
  int2size_t::iterator it = nucsize.find(taxid);
  assert(it!=nucsize.end());
  return it->second;
}


size_t getsum_of_subtree(int2size_t &retmap, int2intvec &child, int taxid) {
  size_t ret = 0;
  int2size_t::iterator it = retmap.find(taxid);
  if (it != retmap.end())//if we are at leafnode or we have computed the sum of subtree
    ret += it->second;
  //  fprintf(stderr,"retexists\t%d\t%lu\n",taxid,ret);
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
    //maybe below doesnt work
    //  fprintf(stderr,"retdone\t%d\t%lu\n",taxid,ret);
    if(it==retmap.end())
      retmap[taxid] = ret;
    else
      it->second = ret;
    return ret;
}

size_t how_many_subnodes(int2intvec &child,int taxid){
  assert(taxid!=-1);

  int2intvec::iterator it = child.find(taxid);
  if(it==child.end())
    return 1;
  std::vector<int> &avec = it->second;
  if(avec.size()==0){ //I dont think this can happen{
    fprintf(stderr,"avec.size() is never zero\n");
    exit(1);
    return 1;
  }
  
  size_t nsum =1;
  for(int i=0;i<avec.size();i++)
    if(avec[i]!=-1){
      nsum += how_many_subnodes(child,avec[i]);
    }
  return nsum;
}


void print_a_subtree(FILE *fp, int2intvec &child,int taxid){
  assert(taxid!=-1);
  fprintf(fp,"%d\n",taxid);
 
  int2intvec::iterator it = child.find(taxid);
  if(it!=child.end()){
    std::vector<int> &avec = it->second;
    if(avec.size()==0)
      fprintf(stderr,"I dont think this can happen \n");
      
    for(int i=0;i<avec.size();i++)
      if(avec[i]!=-1)
	print_a_subtree(fp,child,avec[i]);
  }
}

//this small function will loop over all ranks for all taxids
//if a rank is -1, then it is called 'no rank'
// a 'no rank' can be found anywhere in our tree structure.
//This function will set the the level of a no rank to subspecies level if it is below species. And will set it to kingdom (or higher) level if it is above species level.
//the rationale is that we contract all nodes below species to species
//and discard those above. 
void fix_no_rank(int2int &rank,int2int &parent){
  fprintf(stderr,"Calling fix rank\n");
  int nadjust = 0;
  for(auto &x:rank){
    if(x.second!=-1)
      continue;
    nadjust++;
    int up_has_species  = 0;
    int2int::iterator it = parent.find(x.first);
    assert(it!=parent.end());
    while(it->second!=1){
     
      int2int::iterator itr = rank.find(it->second);
      assert(itr!=rank.end());
      if(itr->second==5){
	up_has_species = 1;
	break;
      }else if(itr->second>5)
	break;
      it = parent.find(it->second);
      assert(it!=parent.end());
    }
    if(up_has_species ==1)
      x.second = 5;
    else
      x.second = 37;
  }
  fprintf(stderr,"Done calling fix rank. We have adjusted: %d nodes propotion is: %f\n",nadjust,(double) nadjust/(double) rank.size());
}
/*
  this function should be called from root node.
  Whenever a node with species level has been found then the remaining subtree will be removed from main tree.
 */

void prune_below_species(int2int &parent, int2intvec &child,int2int &rank,int taxid){
  assert(taxid!=-1);

  int2int::iterator iti = rank.find(taxid);
  assert(iti!=rank.end());
  int myrank = iti->second;
  int2intvec::iterator itv = child.find(taxid);

  if(myrank>5){//rank==5 is species

    //only if we have child nodes
    if(itv!=child.end()){
      std::vector<int> &avec = itv->second;
      assert(avec.size()>0);
      for(int i=0;i<avec.size();i++){
	if(avec[i]==-1)
	  continue;
	prune_below_species(parent,child,rank,avec[i]);
      }
    }
  }else if(myrank==5) {
    //only if we have childnodes
    if(itv!=child.end()){
      std::vector<int> &avec = itv->second;
      assert(avec.size()>0);
	
      for(int i=0;i<avec.size();i++){
	if(avec[i]==-1)
	  continue;

	//first set the node from the parent of the downnode to =-1
	int2int::iterator it2 = parent.find(avec[i]);
	if(it2==parent.end()){
	  fprintf(stderr,"Problem finding taxid: %d from taxid: %d\n",avec[i],taxid);
	  exit(1);
	}
	it2->second = -1;
	//then set the downnode of the parent to =-1;
	avec[i] = -1;
      }
    }
  }else if(myrank<5){
    fprintf(stderr,"I dont think this can happen myrank<5: taxid: %d \n",taxid);
  }else{
    fprintf(stderr,"Never here\n");
  }
  

}


void print_leafs(FILE *fp, int2intvec &child,int taxid,int2size_t &hasgenome){
  assert(taxid!=-1);//shouldnt happen
 
  int2intvec::iterator it = child.find(taxid);
  if(it==child.end()){
    int2size_t::iterator its = hasgenome.find(taxid);
    if(its!=hasgenome.end())//only print if leaf with data
      fprintf(fp,"%d\t%lu\n",taxid,its->second);
  }else{
    std::vector<int> &avec = it->second;
    if(avec.size()==0)
      fprintf(stderr,"I dont think this can happen \n");
    
    for(int i=0;i<avec.size();i++)
      if(avec[i]!=-1)
	print_leafs(fp,child,avec[i],hasgenome);
    
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
	    if(kstr->s[0]=='#')
	      continue;
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

//this function assumes that there is only data at leafnodes
//taxid is topnode of subtree
//child is map of child nodes
//num_subgenomes is map of number of genomes is subtree
//nrep is the number of representative genomes we request
//ret will contain an array of genomes that has been selected as being representative
void getnrep(int taxid,int2intvec &child,int2size_t &num_subgenomes,int nrep,std::vector<int> &ret){
  if(nrep==0)
    return;
  int2size_t::iterator it = num_subgenomes.find(taxid);
  assert(it!=num_subgenomes.end());//validate that the data structure exists
  if(it->second==0){
    fprintf(stderr,"We should never be in the situation where we dont have any genomes available in clade\n");
    return;
  }
  size_t numG = it->second;
  
  //"We should never be in the situation where we request more genomes than what we have available\n");
  assert(nrep<=numG);
  

 
  int2intvec::iterator it2 = child.find(taxid);
  if(it2==child.end()){//this is leaf
    //  fprintf(stderr,"\t\t-> Adding taxid: %d, done with this call\n",taxid);
    assert(numG==1);
    ret.push_back(taxid);//this adds the taxid to the results
    return;
  }
  
  // we have child nodes
  std::vector<int> &avec = it2->second;
  int nchild = 0;
  for (unsigned i = 0; i < avec.size(); i++)
    if(avec[i]!=-1&&getval(num_subgenomes,avec[i])>0)
      nchild++;

  //nchild is the number of childs that has data

  assert(nchild>0);

  int at =0;
  int ntaken = 0;
#if 0
  for (unsigned i = 0; i < avec.size(); i++) {
    if(avec[i]!=-1&&getval(num_subgenomes,avec[i])>0)
      fprintf(stderr,"INFO: taxid:%d avec:%d subG:%lu numG:%lu child:%d\n",taxid,avec[i],getval(num_subgenomes,avec[i]),numG,nchild);
  }
#endif
  for (unsigned i = 0; i < avec.size(); i++) {
    if(avec[i]!=-1&&getval(num_subgenomes,avec[i])>0){
      
      int target = ceil((double)(nrep-ntaken)/((double)nchild-at));
      //validate that we are not requesting more genomes than what we have available
      if(target>getval(num_subgenomes,avec[i]))
	target = getval(num_subgenomes,avec[i]);
      
      numG -= getval(num_subgenomes,avec[i]);
          
      while(numG+target<nrep-ntaken){
	target++;
      }

      //      fprintf(stderr,"%d) taxid: %d target: %d ngenom:%lu NEWNODE: %D\n",at,taxid,target,numG,avec[i]);
      at++;
      ntaken += target;
      getnrep(avec[i],child,num_subgenomes,target,ret);;
    }
  }
 
}

int *splitdb(int nk,int2int &up,int2intvec &down,int2size_t &nucsize){
  size_t total_sum_bp = nucsize.find(1)->second;
  size_t total_sum_nodes = how_many_subnodes(down,1);
  double target = total_sum_bp/nk;
  fprintf(stderr,"\t->[%s] total sum of tree: %lu target for subtrees: %f total_nodes: %lu\n",__FUNCTION__,total_sum_bp,target,total_sum_nodes);
  int *trees = new int[nk];

  for(int k=0;k<nk-1;k++){//loop over subtrees
    fprintf(stderr,"\t->------------ k:%d nk:%d howmany_nodes: %lu totalbp: %lu\n",k,nk,how_many_subnodes(down,1),getval(nucsize,1));
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

char2int getrank(){
  char2int ret;
  
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


size_t get_total_length2(faidx_t *myfai){
  size_t total = 0;
  for(int i=0;i<faidx_nseq(myfai);i++){
    const char *name = faidx_iseq(myfai,i);
    total += faidx_seq_len64(myfai,name);
  }
  return total;
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
    gzclose(fp2);
  }

  gzclose(fpout);
  gzclose(fp);
  for(auto x: ass2tax) free(x.first);

  return 0;
}

size_t read_taxid_bp(char *fname,int2size_t &ret){
  assert(fname!=NULL);
  BGZF *fp = NULL;
  assert(((fp=bgzf_open(fname,"rb")))!=NULL);

  kstring_t kstr;
  kstr.l=kstr.m=0;
  kstr.s=NULL;
  int2size_t::iterator it;
  while (bgzf_getline(fp, '\n', &kstr)) {
    if(kstr.l==0)
      break;
    //      fprintf(stderr,"line: %s\n",kstr.s);
    int taxid = atoi(strtok(kstr.s,"\t\n "));
    if(taxid<1)
      continue;
    size_t bp = atol(strtok(NULL,"\t\n "));
    it = ret.find(taxid);
    if(it==ret.end())
      ret[taxid] = bp;
    else
      it->second = it->second +  bp;
    kstr.l = 0;
  }
  bgzf_close(fp);
  free(kstr.s);
  return ret.size();
}


int main_filter(char *fname,int2int &parent,int2intvec &child,int2int &rank){
  fprintf(stderr,"[%s]\n",__FUNCTION__);
  int2size_t ret;
  read_taxid_bp(fname,ret);
  for(auto x: ret){
    int2int::iterator it = parent.find(x.first);
    int2intvec::iterator itv = child.find(x.first);
    int2int::iterator iti = rank.find(x.first);
    

    //    assert(itv!=child.end());
    if(it==parent.end()){
      fprintf(stderr,"ERR_noparent\t%d\t%lu\n",x.first,x.second);
      continue;
    }

    assert(iti!=rank.end());//<- this should be done after checking for parent 
   
    if(iti->second > 5){
      fprintf(stderr,"ERR_rank_too_high\t%d\t%lu\n",x.first,x.second);
      continue;
    }
    if(iti->second == 5){
      fprintf(stdout,"%d\t%u\n",x.first,x.second);
      continue;
    }
    //we have now checked those with correct rank, and discarded those with too high.
    //for those that are subspecies or similar.
    int newtaxid = x.first;
    fprintf(stderr,"ERR_toolow_rank_will_loopup:original_%d->",newtaxid);
    while(iti->second<5){
      if(iti->second>5){
	fprintf(stderr,"Big problem rank of taxid is resolved at too high level\n");
      }
      it = parent.find(newtaxid);
      assert(it!=parent.end());
      newtaxid = it->second;
      iti = rank.find(newtaxid);
      assert(iti!=rank.end());
      fprintf(stderr,"%d->",newtaxid);
    }
    //validate that the next node is indeed above species level, (species group)
    int2int::iterator iti_extra = parent.find(newtaxid);
    assert(iti_extra!=parent.end());
    int2int::iterator itv_extra = rank.find(iti_extra->second);
    assert(itv_extra!=rank.end());
    if(itv_extra->second<=5){
      newtaxid = iti_extra->second;
      iti = rank.find(newtaxid);
    }
    //do sanity check again. This can be removed
    iti_extra = parent.find(newtaxid);
    assert(iti_extra!=parent.end());
    itv_extra = rank.find(iti_extra->second);
    assert(itv_extra!=rank.end());
    assert(itv_extra->second>5);
    //delete above 5 lines. 
      
    fprintf(stderr,"Done Looping, final rank: %d\n",iti->second);
    fprintf(stdout,"%d\t%u\n",newtaxid,x.second);
    
  }
  return 0;
}




int main(int argc,char **argv){
  char *node_file = strdup("nodes_v6.dmp.gz");
  char *meta_file = "/projects/lundbeck/scratch/taxDB/v6/metadata/taxdb-genome_stats-broad-v6.tsv.gz";
  const char *outname = "outname";
  char *acc2taxid_flist = "acc2taxid_flist";
  int how_many_chunks = 8;
  int makefai = 0;
  char *wgs_fname = NULL;
  char *seqs_fname = NULL;
  int filter = 0;
  int nrep = 10;
  for(int at = 1;at<argc;at++){
    if(strcasecmp(argv[at],"-node_file")==0){
      free(node_file);
      node_file = strdup(argv[at+1]);
    }
    else if(strcasecmp(argv[at],"-meta_file")==0)
      meta_file = strdup(argv[at+1]);
    else if(strcasecmp(argv[at],"-acc2taxid_flist")==0)
      acc2taxid_flist = strdup(argv[at+1]);
    else if(strcasecmp(argv[at],"-nchunks")==0)
      how_many_chunks = atoi(argv[at+1]);
    else if(strcasecmp(argv[at],"-nrep")==0)
      nrep = atoi(argv[at+1]);
    else if(strcasecmp(argv[at],"-outname")==0)
      outname = strdup(argv[at+1]);
    else if(strcasecmp(argv[at],"makefai")==0){
      makefai = atoi(argv[at+1]);
    }
    else if(strcasecmp(argv[at],"filter")==0){
       filter = atoi(argv[at+1]);
    }
    else if(strcasecmp(argv[at],"-wgs")==0){
      wgs_fname = strdup(argv[at+1]);
    }
    else if(strcasecmp(argv[at],"-seqs")==0){
      seqs_fname = strdup(argv[at+1]);
    }
    at++;;
      
  }
  fprintf(stderr,"\t 1) ./program -node_file filename.txt -nchunks integer -nrep %d -wgs %s -seqs %s\n\t 2) ./program makefai 1 -meta_file filenames.txt -acc2taxid_flist file.list\n\t 3) ./program filter 1 -meta_file taxid_bp.txt -node_file filename.txt \n\t-> -node_file: \'%s\'\n\t-> -meta_file: \'%s\'\n\t-> -nchunks: %d\n\t-> -outname: %s\n\t-> -acc2taxid_flist: %s\n\t-> makefai: %d\n\t-> filter: %d\n",nrep,wgs_fname,seqs_fname,node_file,meta_file,how_many_chunks,outname,acc2taxid_flist,makefai,filter);
  if(makefai){
    return main_fai_extension(meta_file,outname,acc2taxid_flist,makefai);
  }
  
  int2int taxid_rank;
  int2int taxid_parent;
  int2intvec taxid_childs;

  parse_nodes(node_file,taxid_rank,taxid_parent,taxid_childs,1);
  fix_no_rank(taxid_rank,taxid_parent);
  //  print_node(stderr,taxid_parent,taxid_childs,1);
  fprintf(stderr,"\t-> howmany nodes: %lu taxid_parents.size(): %lu, these values should be identical\n",how_many_subnodes(taxid_childs,1),taxid_parent.size());
  assert(how_many_subnodes(taxid_childs,1)==taxid_parent.size());

  if(filter){
    return main_filter(meta_file,taxid_parent,taxid_childs,taxid_rank);
  }
  //  print_node(stderr,taxid_parent,taxid_childs,46014);

  prune_below_species(taxid_parent,taxid_childs,taxid_rank,1);
  
#if 0
  for(int2int::iterator it=taxid_parent.begin();it!=taxid_parent.end();it++)
    fprintf(stderr,"\t%d) -> %d\n",it->first,it->second);
#endif
  //  char2int ass2tax=acc2taxid(acc2taxid_flist);

 
  int2size_t wgs_map,tmp_map,total_map;
  //read in wgsfiles
  read_taxid_bp(wgs_fname,wgs_map);
  total_map=tmp_map = wgs_map;
  //compute aggregated sum
  getsum_of_subtree(wgs_map, taxid_childs, 1);
  fprintf(stderr,"\t-> Number of basepairs from wgs: %lu\n",getval(wgs_map,1));

  
  //read in seqfiles
  if(seqs_fname!=NULL)
    read_taxid_bp(seqs_fname,total_map);
  getsum_of_subtree(total_map, taxid_childs, 1);
  fprintf(stderr,"\t-> Number of basepairs from wgs+seqs: %lu\n",getval(total_map,1));
  //  read_taxid_bp(seqs_fname,total_map);

  // int2size_t taxid_genome_size =parse_meta(meta_file);
  int2size_t::iterator it = total_map.begin();
  fprintf(stderr,"\t-> total presize: %lu key: %d val:%lu\n",total_map.size(),it->first,it->second);

  fprintf(stderr,"\t-> postsize: %lu ratio: %f \n",total_map.size(),(float)total_map.size()/wgs_map.size());
  //  exit(0);
#if 0
  for(int2size_t::iterator it=total_map.begin();it!=total_map.end();it++)
      if(it->second!=0)
      fprintf(stderr,"\t%d) -> %lu\n",it->first,it->second);
  return 0;
#endif
  
  int *subtrees = splitdb(how_many_chunks,taxid_parent,taxid_childs,total_map);

  //subtrees contains the taxids of the "root" of each detached subtree
  for(int i=0;i<how_many_chunks;i++){
    //    fprintf(stderr,"chunk:%d taxid: %d \thow_many:%lu\tget_val:%lu\n",i,subtrees[i],how_many_subnodes(taxid_childs,subtrees[i]),getval(total_map,subtrees[i]));
    //fprintf(stderr,"sub[%d]: %d\n",i,subtrees[i]);
    char onam[1024];
    snprintf(onam,1024,"%s_cluster.%d",outname,i);
    FILE *fp = fopen(onam,"wb");
    fprintf(stderr,"\t-> Writing file: %s\n",onam);
    print_leafs(fp,taxid_childs,subtrees[i],tmp_map);
    fclose(fp);
  }
  
  wgs_map = tmp_map;
  //wgs_map contains the size of the different genomes.
  //we will not change this so it is simply an indicator of whether or not we have data
  for(int2size_t::iterator it=wgs_map.begin();it!=wgs_map.end();it++)
    it->second = 1;

  //for each subtree calculate the sum of genomes
  for(int i=0;i<how_many_chunks;i++)
    getsum_of_subtree(wgs_map,taxid_childs,subtrees[i]);
  
  for(int i=0;i<how_many_chunks;i++){
    char onam[1024];
    snprintf(onam,1024,"%s_representative.%d",outname,i);
    FILE *fp = fopen(onam,"wb");
    fprintf(stderr,"\t-> Writing file: %s\n",onam);
    std::vector<int> genomes;
    getnrep(subtrees[i],taxid_childs,wgs_map,nrep,genomes);
    for(int j=0;j<genomes.size();j++){
      int2size_t::iterator its = tmp_map.find(genomes[j]);
      assert(its!=tmp_map.end());
      fprintf(fp,"%d\t%lu\n",genomes[j],its->second);
    }
    fclose(fp);
  }

  free(node_file);
  delete [] subtrees;
  free(wgs_fname);
  free(seqs_fname);
  return 0;
}
