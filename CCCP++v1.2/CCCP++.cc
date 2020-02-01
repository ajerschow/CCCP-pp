#include<stdio.h>
#include<math.h>
//#include "qsim.h"
#include <FlexLexer.h>
#include <vector>
#include <list>
#include <iostream>
#include "ParsedVars.h"
#include "gcd.h"
#include <algorithm>  // Include algorithms
#include <string>

#ifdef BRUKEROUTPUT
#include <sstream>
#endif

#include <fstream>

#define VERSION "1.2 (2005)"

using namespace std;
// contains max. number in ph
extern int Nph;
vector<float> PCsig;
//using namespace std;


// cog reading and verifying apparently correct
// newattn contains the deviation coh. order
// next: go ahead with COGwheel algo



// tested, lexer works with qsim parameters, reads into arrays.
// reads list of vectors correctly
// next: read in real PC parameters, convert PC, gcd, scm
// for PC, can't one calculate correlation or anticorrelation between
// different phases? multiple phases? linear comb?
// what's the best way to read in attn (for each nucleus)?

// Symmetry: negative pathways exact same thing up to receiver phase
// higher symmetry possible?
// +- for each step (each pulse) -> maybe too complicated, more compl. than adv.
// calc. phase up to receiver, there choose between positive and negative ones.
// Add to list.




int safemod(int a, int N){
  if (a>=0)
    return a%N;
  else
    if (a==(a/N)*N)
      return 0;
    else
      return a+(-a/N+1)*N;
}


bool SRTpw(int a,int b){
  return PCsig[a]<PCsig[b];
}

template <class T>
void reset(vector<T>& v,T val){
  for (int i=0;i<v.size();i++){
    v[i]=val;
  }
}

inline float abs(float f)
{
  return f >= 0.0 ? f: -f;
}
 
template <class T>
ostream& operator<< (ostream& s, const vector<vector<T> >& A) {
  int i,j;
  for (i=0;i<A.size();i++){
    for (j=0;j<A[i].size();j++){
      cout << A[i][j] << " ";
    }
    cout<<endl;
  }
  cout<<endl;
}


template <class T>
ostream& operator<< (ostream& s, const list<vector<vector<T> > >& A);

template <class T>
ostream& operator<< (ostream& s, const list<vector<vector<T> > >& A) {
  typename list<vector<vector<T> > >::iterator lattn;
  for (lattn=A.begin();lattn!=A.end();lattn++){
    cout << (*lattn);
  }
}


template <class T>
ostream& operator<< (ostream& s, list<list<vector<T> > >& A) {
  typename list<list<vector<T> > >::iterator lattn;
  for (lattn=A.begin();lattn!=A.end();lattn++){
    cout << (*lattn);
  }
}

// includes sorted output
void outputpaths(vector<list<vector<int> > > & paths,vector<int> pidx,vector<float> PCsig,vector<float> PCphase){
  int nucind;
  list<vector<int> >::iterator ip;
  for (int i=0;i<paths.size();i++){
    cout<<"#"<<i<<":\n";
    nucind=1;
    for(ip=paths[pidx[i]].begin();ip!=paths[pidx[i]].end();ip++){
      cout<<"  nuc "<<nucind<<": ";
      for (int ii=0;ii<(*ip).size();ii++){
	cout<<((*ip)[ii])<<" ";
      }
      cout<<endl;
      nucind++;
    }
    cout << "  signal: "<<PCsig[pidx[i]]<<" phase: "<<PCphase[pidx[i]]<<endl;
  }
}

// this one adds COGpw back in
void outputpaths(vector<list<vector<int> > > & paths, list<vector<int> > COGpw){
  int nucind;
  list<vector<int> >::iterator ip,cp;
  for (int i=0;i<paths.size();i++){
    cout<<"#"<<i<<":\n";
    cp=COGpw.begin();
    nucind=1;
    for(ip=paths[i].begin();ip!=paths[i].end();ip++){
      cout<<"nuc "<<nucind<<": ";
      for (int ii=0;ii<(*ip).size();ii++){
	cout<<((*ip)[ii]+(*cp)[ii])<<" ";
      }
      cout<<endl;
      cp++;
      nucind++;
    }
  }
}

template <class T>
ostream& operator<< (ostream& s, list<vector<T> >& A) {
  int i,j;
  typename list<vector<T> >::iterator il;
  for (il=A.begin();il!=A.end();il++){
    for (j=0;j<(*il).size();j++){
      cout << (*il)[j] << " ";
    }
    cout<<endl;
  }
  cout<<endl;
}

ostream& operator<< (ostream& s,  list<int>& A) {
  int i,j;
  list<int>::iterator il=A.begin();
  for (il=A.begin();il!=A.end();il++){
    cout << (*il) <<" " " ";
  }
  cout<<endl;
}

template <class T>
ostream& operator<< (ostream& s, const vector<T>& A) {
  int i,j;
  for (i=0;i<A.size();i++)
    cout << A[i] << " ";
  cout<<endl;
}

template <class T>
void operator/= (vector<T>& A, T a) {
  int i;
  for (i=0;i<A.size();i++)
    A[i]/=a;
}

template <class T>
vector<T> operator* (const vector<T>& A, T a) {
  int i;
  vector<T> res(A.size());
  for (i=0;i<A.size();i++)
    res[i]=A[i]*a;
  return res;
}

template <class T>
vector<T> operator- (const vector<T>& A) {
  int i;
  vector<T> res(A.size());
  for (i=0;i<A.size();i++)
    res[i]=-A[i];
  return res;
}

template <class T>
int allequal(const vector<T>& A) 
{
  int i;  T el=A[0];
  if (A.size()==1)
    return 1;
  else
    {
      i=1;
      while (A[i]==A[i-1] && i<A.size())
	i++;
      if (i==A.size())
	return 1;
      else
	return 0;
    }
}


int condition_phase(list<vector<int> >& ph, vector<vector<int> >& newphase, vector<int>& receiver_ph, vector<int>& base, int & fullbase, vector<int>& phase_allequal)
{
  int i,j,ii,dummy;
  list<vector<int> >::iterator phasel;

  // better to eliminate lines with constant phases, eliminate also in attn
  // need to calc. everything mod base first.
  // keep track of removed lines
  // subtract const. phase, so can remove first step?
  //     don't forget to add it back in when calc. attenuation
  // implement exponential better
  // ORGANIZE YOURSELF

  // think about partitions, inhomog. and homog. solution
  // do partitions help?

  /////////////////////////////////////////////////
  ////////////////////////////////////
  // CAREFUL ABOUT GCD ROUTINES     //
  // GCD SHOULD INCLUDE BASE??????  YEAH  //
  ////////////////////////////////////

  if (base.size()==1)
    for (ii=0;ii<(ph.size());ii++)
      base.push_back(base[0]);
  else if (base.size()!=ph.size()){
    cerr << "base vector does not match the size of the phase cycle\nAborting...\n";
    exit(1);
  }

#ifdef VERBOSEXXX
  cout << "Base " << base << endl;
#endif

  phasel=ph.begin();
  for (i=0;i<ph.size();i++)
    {
      //issue warning if size not divider of Nph
      if (Nph % (*phasel).size() !=0)
	cerr << "Warning, phase " << i+1 << " is has a length that does not divide " << Nph << endl; 

      // make all phases mod the appropriate base
      for (j=0;j<(*phasel).size();j++)
	(*phasel)[j]=safemod((*phasel)[j],base[i]);

      // modify base to reflect minimal base:
      // properly deal with zero bases, no change in phase
      dummy=gcd(*phasel);
      dummy=gcd(dummy,base[i]);
#ifdef VERBOSEXXX
      cout <<  (*phasel) << "gcd " << dummy << endl;
#endif
      if (dummy!=0)
	{
	  base[i]=base[i]/dummy;
	  (*phasel)/=dummy;
	}
      else
	base[i]=1;
      phasel++;
    }

#ifdef VERBOSEXXX
  cout << "Base " << base << endl;
#endif
  // find common base
  fullbase=lcm(base);


  // GET RID OF EFFICIENCY-MADNESS....
  // NO OPTIM. OF PC
  // CAN BE DONE LATER USING BASE PROBING, ETC.
  // after this, simply remove all base=1 phases, but keep track
  // and subtract first phase
  // put into vector<vector>, store first_phase, store base[i] for future, 
  // store phase_allequal..... 
  
  // maybe better first remove base=1 phases, so that fullbase may be smaller!
  // on the other hand fullbase does not really lead to loss in efficiency

  //cout << "fullbase "<<fullbase<<endl;
  vector<int> dummyvec;
  phasel=ph.begin();

  for (i=0;i<ph.size();i++)
    {
      // record which phases
      phase_allequal.push_back(allequal(*phasel));
	    
      dummyvec.clear();
      for (j=0;j<Nph;j++){
	//cout << (((*phasel)[j%(*phasel).size()] -dummy)%base[i]) *fullbase/base[i]<< " ";
	//dummyvec.push_back( ( ((*phasel)[j%(*phasel).size()]-dummy)%base[i]) *fullbase/base[i]);
	//dummyvec.push_back( safemod((*phasel)[j%(*phasel).size()]-dummy,base[i])*fullbase/base[i]);
	dummyvec.push_back( safemod((*phasel)[j%(*phasel).size()],base[i])*fullbase/base[i]);
      }
      //cout << endl;
      newphase.push_back(dummyvec);
      phasel++;
    }


  //// pointer phasel now points to last element, i is last index
  //// no minus since the whole thing is turning the other way -Delta p * phi - phi_r
  //// so this hardly has any meaning
  // // drop off newphase;
  receiver_ph=newphase[newphase.size()-1];
  newphase.pop_back();
  cout << "###########################################################" << endl;
  cout << "# Calculate selected coherence pathways by phase cycle"<< endl;
  cout << "###########################################################" << endl;

  cout << "Phase cycle steps: " << Nph << endl;
  cout << "adjusted phases: " << endl << newphase;
  cout << "receiver: " << receiver_ph<<endl;
  cout << "base: "<< base<<endl;

}

void increment_vec(vector<int>& v, const vector<int> &max)
{
  int inc_idx=0;
 
  do {
    if((++v[inc_idx])>=max[inc_idx]){
      v[inc_idx++]=0;
      v[inc_idx]++;
    }
  } while (v[inc_idx]>=max[inc_idx]);
}


void calc_Deltap(float & abssig, float & sigphase, vector<int> & diffpw, vector<vector<int> >& newphase,vector<int>& receiver_ph, vector<float> & sin_st, vector<float>& cos_st,float & cutoffSQ, int & fullbase)
{
  // if necessary, create transpose of newphase mtx first
  // would be interesting to see how much this gives
  // possib. speedup: store sin/cos, make angle power of 2, possible?
  // so modulo becomes AND operation

  float sig_c,sig_s;

  // receiver
  vector<int> dp_ph=receiver_ph;
  // cout << "rec: " << receiver_ph;
  for (int idx2=0;idx2<newphase[0].size();idx2++){
    for (int idx1=0;idx1<newphase.size();idx1++){
      dp_ph[idx2]+=diffpw[idx1]*newphase[idx1][idx2];
    }
    dp_ph[idx2]=safemod(dp_ph[idx2],fullbase);
  }

  //cout << "dp_ph " << dp_ph;

  // calc. signal
  sig_c=0.0;sig_s=0.0;
  for (int i=0;i<Nph;i++){
    //cout << "(" << sin_st[dp_ph[i]] << " " << cos_st[dp_ph[i]] << ") ";  
    sig_s+=sin_st[dp_ph[i]];
    sig_c+=cos_st[dp_ph[i]];
  }

  abssig=sig_s*sig_s+sig_c*sig_c;
  if (abssig>cutoffSQ){
    sigphase=atan2(sig_s,sig_c);
    abssig=sqrt(abssig);

  }
  else {
    abssig=0.0;
    sigphase=0.0;
  }


}


list<vector<vector<int> > >  convert_attn(list<vector<vector<float> > >& attn, list<vector<int> > & newattn_max, vector<int> & which_nucl, vector<int> & nucl_pw_idx)
{
  // this one stores the result, which is a list for each nucleus
  // and a vector of vectors for each pulse, returning also the maximum number of coherences
  // for each pulse.
  // test if which_nucl compatible
  // which_nucl should be [0,maxnucl-1];
  // should now be done in main()
  // 
  list<vector<vector<int> > > res;
  int i,j,startind,maxval;
  // stores temporarily the coh. orders for one pulse
  vector<int> dummyvec,dummymaxvec;
  // stores temporarily all coh. order vectors for one nucleus
  vector<vector<int> > dummyvecvec;



  vector<int> countpulses(attn.size(),0);
  nucl_pw_idx.clear();
  for (i=0;i<which_nucl.size();i++){
    if (which_nucl[i]<0 || which_nucl[i]>(attn.size()-1)){
      cerr << "which_nucl contains an out of range nucleus index "<< which_nucl[i]+1 << " at pulse " << i <<"!\nAborting....\n";
      exit(1);
    }
    nucl_pw_idx.push_back(countpulses[which_nucl[i]]);
    ++countpulses[which_nucl[i]];
  }

  newattn_max.clear();

  list<vector<vector<float> > >::iterator nucit;
  list<vector<int> >::iterator COGpwit=COGpw.begin();

  int  nucind=0;

  // here: for COG transform newattn such that only the difference 
  // between the desired pw and the others is cycled through
  // at the same time can check if COGpw has the right size
  if ((COGmax_pw!=NOTASSIGNED && COGmax_N!=NOTASSIGNED)){
    if (attn.size()!=COGpw.size()) {
      cerr<< "'COGpw' is not for the same number of nuclei as 'attn', i.e. " <<attn.size()<<endl;
      exit(1);
    }
  }

  // considers attn made of length appropriate for given nucleus, is this good?
  for (nucit=attn.begin(); nucit!=attn.end();nucit++){
    if ((*nucit)[0].size()!=(countpulses[nucind]+1)){
      cerr << "attn["<<nucind+1<<"] does not match the number of\npulses indicated in which_nucl!\nShould be "<<countpulses[nucind]+1<<"\nAborting...\n";
      exit(1);
    }

    // check size of COGpw
    if (COGmax_pw!=NOTASSIGNED && COGmax_N!=NOTASSIGNED){
      if ((*COGpwit).size()!=(countpulses[nucind]+1)){
	cerr << "COGpw["<<nucind+1<<"] does not match the number of\npulses indicated in which_nucl!\nShould be "<<countpulses[nucind]+1<<"\nAborting...\n";
	exit(1);
      }
    }

    dummyvecvec.clear();
    dummymaxvec.clear();
    startind=((*nucit).size()-1)/2;

    for (j=0;j<(*nucit)[0].size();j++){
      dummyvec.clear();
      maxval=0;
      for (i=0;i<(*nucit).size();i++){
	if ((*nucit)[i][j]!=0.0){
	  // add in newattn the difference pw betw. actual and desired
	  if (COGmax_pw==NOTASSIGNED || COGmax_N==NOTASSIGNED){
	    dummyvec.push_back(startind-i);}
	  else {
	    dummyvec.push_back(startind-i-(*COGpwit)[j]);
	  }
	  maxval++;
	}
      }
      if (maxval==0){
	cerr << "attn["<< nucind+1 <<"] does not allow any pathways, at pulse "<< j+1 << "!\n Aborting..." << endl;
	exit(1);
      }
      dummymaxvec.push_back(maxval);
      dummyvecvec.push_back(dummyvec);
    }
    res.push_back(dummyvecvec);
    newattn_max.push_back(dummymaxvec);
    nucind++;
    COGpwit++;
  }
  return res;
}


/////////////////////////////////////////////
///      REWRITE as LIST<VECTOR,INT> WITH MAX
///      CHECKING, SO IT'S FASTE/////////////////////////////////////////////


void initialize_pw_idx(list<vector<int> > & pw_idx,  list<vector<vector<int> > > & newattn)
{
   pw_idx.clear();
   list<vector<vector<int> > >::iterator newattn_it;
   //   list<vector<int> >::iterator list_it;
   vector<int> dummyvec;
   
   for (newattn_it=newattn.begin();newattn_it!=newattn.end();newattn_it++) {
     dummyvec.clear();
     for (int i=0;i<(*newattn_it).size();i++)
       dummyvec.push_back(0);
     pw_idx.push_back(dummyvec);
   }
}
    
int increment_pw_idx(list<vector<int> > & pw_idx,  list<vector<vector<int> > >& newattn, list<vector<int> > & newattn_max)
{
  list<vector<int> >::iterator pw_idx_it, newattn_max_it;
  list<vector<vector<int> > >::iterator newattn_it;

  int i;
  pw_idx_it=pw_idx.begin();
  newattn_max_it=newattn_max.begin();
  newattn_it=newattn.begin();
  int currpos=0;

  while( (pw_idx_it!=pw_idx.end()) && (++(*pw_idx_it)[currpos]==(*newattn_max_it)[currpos]) ){
    (*pw_idx_it)[currpos++]=0;
    if (currpos==(*pw_idx_it).size()) {
      currpos=0;
      pw_idx_it++;
      //cout <<"pw_idx"<< pw_idx;
      newattn_max_it++;
      //cout << "newattn_max_it "<<(*newattn_max_it)<<endl;
      newattn_it++;
      //cout << "newattn_it"<<(*newattn_it)<<endl;
    }
  } 
  return (pw_idx_it!=pw_idx.end());
}
 

// obtains the current pw (in numbers) from the list of indices into newattn
// does this account for newattn not being of equal length?
void make_pw(list<vector<int> >& currpw, list<vector<int> > & pw_idx, list<vector<vector<int> > > & newattn)
{
  list<vector<int> >::iterator pw_idx_it, currpw_it;
  list<vector<vector<int> > >::iterator newattn_it;

  int i;
  pw_idx_it=pw_idx.begin();
  currpw_it=currpw.begin();
  newattn_it=newattn.begin();

  for (int i=0;i<pw_idx.size();i++){
    for (int j=0;j<(*pw_idx_it).size();j++){
      (*currpw_it)[j]=(*newattn_it)[j][(*pw_idx_it)[j]];
    }
  
    pw_idx_it++;
    currpw_it++;
    newattn_it++;
  }
}

void make_diffpw(vector<int>& diffpw, list<vector<int> >& currpw, vector< list<vector<int> >::iterator >& index_to_currpw, vector<int>& which_nucl, vector<int>& phase_allequal, vector<int>& nucl_pw_idx)
{
  list<vector<int> >::iterator currpw_it;
  vector<int>::iterator phase_allequal_it;
  
  int ind=0,pulsind=0,dummy,nucnum,dumidx;

  for(pulsind=0;pulsind<which_nucl.size();pulsind++){
    nucnum=which_nucl[pulsind];
    dumidx=nucl_pw_idx[pulsind];
    //cout << "["<<pulsind<<","<<nucnum<<","<<dumidx<<"] ";
    diffpw[pulsind]=(*(index_to_currpw[nucnum]))[dumidx+1]-(*(index_to_currpw[nucnum]))[dumidx];
  }
  //cout<<endl;
}

int calc_cogwheel_receiver_phase(vector<int> which_nucl,vector<int > COGwdg, list<vector<int> > COGpw, int fullbase){
  // calc receiver phase (need to have index for each nucleus)
  // here, need to calc. receiver phase
  // evtl. output of all phases?
  list<vector<int> >::iterator COGpwit;
  int nucind=0,sum=0,nucpulsind=0;

  for (COGpwit=COGpw.begin();COGpwit!=COGpw.end();COGpwit++){
    nucpulsind=0;
    for (int pulsind=0;pulsind<which_nucl.size();pulsind++){
      //careful, which_nucl should be now starting at zero

      if (which_nucl[pulsind]==nucind){
	// cout << "["<<pulsind<<";---------";
	sum-=((*COGpwit)[nucpulsind+1]-(*COGpwit)[nucpulsind])*COGwdg[pulsind];
	nucpulsind++;
      }
    }
    nucind++;
  }
  return safemod(sum,fullbase);
}

main(int argc, char** argv){
  int   i,j,maxbase;
  int fullbase;
  vector<int> phase_allequal;
  vector<vector<int> > newphase;
  vector<int> receiver_ph;

  // consider symmetry betw. posit. and negat. paths

  if (argc>1 && !strcmp(argv[1],"-v")) {
    cout << "CCCP++, version " << VERSION << ", by Alexej.Jerschow@nyu.edu"<<endl;
    exit(0);
  }
    


  FlexLexer* lexer = new yyFlexLexer;
  while(lexer->yylex() != 0);
  
  //cout << "Lexer successful \n";

  float cutoffSQ=cutoff*cutoff;

  // shift index to nuclei
  for (int i=0;i<which_nucl.size();i++)
    --which_nucl[i];

  // maxjump is max. coh. jump
  if (maxjump==0 || maxjump==NOTASSIGNED)
    maxjump=attn.size();
  if (maxcoh==0 || maxcoh==NOTASSIGNED)
    maxcoh=(attn.size()-1)/2;

  if (cutoff<1e-8)
    cutoff=1e-1;

  // convert phase cycle appropriately,
  // but only if COG not on
  // otherwise will create wild errors if ph not given
  if (COGmax_pw==NOTASSIGNED || COGmax_N==NOTASSIGNED){
    if (COGN!=NOTASSIGNED){
      // create newphase based on winding numbers
      if (COGwdg.size()!=which_nucl.size()){
	cerr<<"COGwdg does not have the right number of winding numbers, i.e. "<<which_nucl.size()<< " as given in which_nucl.....\n";
	exit(1);
      }
      
      cout <<"##############################################################\n";
      cout <<"# COGWHEEL creation of phase cycle based on winding numbers\n";
      fullbase=COGN;
      Nph=fullbase;
      cout <<"##############################################################\n";

#ifdef BRUKEROUTPUT
      cout << "Writing phase cycle\n";
#endif

      // make proper filename
#ifdef BRUKEROUTPUT
      ostringstream fname;
      fname <<"pc_"<<fullbase<<"_";
      for (int i=0;i<which_nucl.size();i++){
	fname<<COGwdg[i];
	//fname=fname+itoa(i);
	if (i<(which_nucl.size()-1))
	  fname<<"_";
      }
      fname<<".bruk";
      ofstream pc_out(fname.str().begin(), ios::out);
#endif

      // calc receiver phase
      int sum=calc_cogwheel_receiver_phase(which_nucl,COGwdg,COGpw,fullbase);


#ifdef BRUKEROUTPUT
      pc_out<<"; COG "<<fullbase<<" (";
#endif
      cout<<"; COG "<<fullbase<<" (";
      for (int i=0; i<which_nucl.size();i++){
#ifdef BRUKEROUTPUT
	pc_out<<COGwdg[i]<<" ";
#endif
	cout<<COGwdg[i]<<" ";
      }
#ifdef BRUKEROUTPUT
      pc_out<<"; "<<sum<<")\n";
#endif
      cout<<"; "<<sum<<")\n";
      
      // cout <<"base= "<<fullbase<<endl;
      for (int ii=0;ii<(ph.size());ii++)
	base.push_back(base[0]);
      
      vector<int> dummyvec;
      int dummyint;
      for (int i=0;i<COGwdg.size();i++){
	dummyvec.clear();
#ifdef BRUKEROUTPUT
	pc_out << "ph" << (i+1) <<"= (" << fullbase << ") ";
#endif
	cout << "ph" << (i+1) <<"= (" << fullbase << ") ";
	for (int j=0;j<fullbase;j++){
	  dummyint=safemod(j*COGwdg[i],fullbase);
	  dummyvec.push_back(dummyint);
#ifdef BRUKEROUTPUT
	  pc_out << dummyint<<" ";
#endif
	  cout << dummyint<<" ";
	}
	phase_allequal.push_back(COGwdg[i]);
	newphase.push_back(dummyvec);
#ifdef BRUKEROUTPUT
	pc_out<< endl;
#endif
	cout<< endl;
      }

#ifdef BRUKEROUTPUT
      pc_out <<"; Receiver:\n";
      pc_out <<"ph0= 0\n";
      pc_out <<"ph31= 0\n";
      pc_out <<"ph30= ("<<fullbase<<") ";
#endif

      cout <<"; Receiver:\n";
      cout <<"ph0= 0\n";
      cout <<"ph31= 0\n";
      cout <<"ph30= ("<<fullbase<<") ";


      for(int i=0; i<fullbase; i++){
	dummyint=safemod(i*sum,fullbase);
	receiver_ph.push_back(dummyint);
	phase_allequal.push_back(dummyint==0);
#ifdef BRUKEROUTPUT
	pc_out<<dummyint<<" ";
#endif
	cout<<dummyint<<" ";
      }
#ifdef BRUKEROUTPUT
      pc_out<<endl;
      pc_out.close();
#endif
      cout<<endl;

#ifdef VERBOSEXXX
      cout << "Receiver: ";
      cout << receiver_ph<<endl;

      cout << "############################################################\n";
#endif

    }
    else {
      // CONVERT ph to newphase, with minimal bases, etc, fullbase given.
      condition_phase(ph, newphase,receiver_ph,base,fullbase,phase_allequal);
      if (ph.size()!=(which_nucl.size()+1)){
	cerr<<"Number of phases does not match the number of pulses indicated in which_nucl, \n i.e. "<<which_nucl.size()<<" plus one receiver phase!\nAborting.....\n";
	exit(1);
      }

      //cout<<"PH\n"<< newphase<<endl;
    }
  }
  // end convert phase cycle


  list<vector<int> > newattn_max;
  // this also verifies if which_nucl and attn contain the correct number of nuclei
  vector<int> nucl_pw_idx;

  // this also subtracts COGpw from entries in newattn if COG is on
  list<vector<vector<int> > > newattn=convert_attn(attn,newattn_max,which_nucl,nucl_pw_idx);

  //cout << "attn" <<endl<<attn;
  //cout << "newattn" <<endl<<newattn;
  //cout << "newattn_max" <<endl<<newattn_max;
  //cout << "phase_allequal" <<endl<<phase_allequal;
  //cout << "nucl_pw_idx " << nucl_pw_idx << endl;
  




  list<vector<int> > pw_idx;
  initialize_pw_idx(pw_idx,newattn);

  //cout << "pw_idx" <<endl<<pw_idx;

  list<vector<int> > currpw=pw_idx;
  list<vector<int> >::iterator currpw_it;
  vector<int> diffpw(which_nucl.size());
  vector< list<vector<int> >::iterator > index_to_currpw;


  //vectors to hold results
  vector<list<vector<int> > > paths;
  //PCsig must be defined globally for sorting
  //vector<float> PCsig;
  vector<float> PCphase;

  //create index to pw index for all nuclei
  for(currpw_it=currpw.begin();currpw_it!=currpw.end();currpw_it++){
    index_to_currpw.push_back(currpw_it);
  }

  vector<float> sin_st;
  vector<float> cos_st;

  if (COGmax_pw==NOTASSIGNED || COGmax_N==NOTASSIGNED){
    // cout <<"fullbase "<<fullbase<<endl;
    // prestore sin,cos
    // prestore sin and cos factors
      sin_st.resize(fullbase);
      cos_st.resize(fullbase);
    for (int i=0;i<fullbase;i++){
      sin_st[i]=sin(2*M_PI*i/(float)fullbase);
      cos_st[i]=cos(2*M_PI*i/(float)fullbase);
    }
    //cout << sin_st;
    //cout<< cos_st;
  }

  float sig_c,sig_s,abssig,sigphase;
  int countremaining=0;



  /////////////////////////////////////////////////////////
  if (COGmax_pw==NOTASSIGNED || COGmax_N==NOTASSIGNED){
    ////////////////////////////////
    // REGULAR PW CALCULATION HERE
    ////////////////////////////////
    do {
      make_pw(currpw,pw_idx,newattn);
      //cout << currpw;
      // for cogwheel: can disregard last Delta p if starting and ending coh. ord same as
      // for desired pw.?
      make_diffpw(diffpw,currpw,index_to_currpw,which_nucl,phase_allequal,nucl_pw_idx);
      // calc. phase cycling sig
      calc_Deltap(abssig,sigphase,diffpw,newphase,receiver_ph,sin_st,cos_st,cutoffSQ,fullbase);
      
      if (abssig>0.0){
	// cout << "currpw " <<currpw;
	// cout << "Sig, phase: " << abssig << ", "<< sigphase<<endl;
	
	PCsig.push_back(abssig);
	PCphase.push_back(sigphase);
	paths.push_back(currpw);
	
      }
      // cout << "----------------" << endl;
      
    } while (increment_pw_idx(pw_idx,newattn,newattn_max));

    ////////////////////////////////////////////////////////////////
    // SORTING
    ////////////////////////////////////////////////////////////////
    int Npaths=PCphase.size();
    vector<int> pidx(Npaths);
    for (int i=0;i<Npaths;i++)
      pidx[i]=i;
    
    // SRTpw is a fctn for sort criterion (so that the index is sorted)
    sort(pidx.begin(),pidx.end(),SRTpw);
    ////////////////////////////////////////////////////////////////
    
    cout<<"--------------------------------------------------"<<endl;
    cout<<"| Selected Pathways:                             |"<<endl;
    cout<<"--------------------------------------------------"<<endl;

    outputpaths(paths,pidx,PCsig,PCphase);
  } else {
    /////////////////////////////////////////////////////////////
    // COGWHEEL PART
    // ALTERNATIVELY, COULD DO LOOP TO AUTOMATICALLY CHECK THE PW
    /////////////////////////////////////////////////////////////
    int COGN;
    int dummy,sum;
    int maxpossib=1;
    vector<int> nu(which_nucl.size());

    // COGwdg can specify which nuclei to avoid cycling
    if (COGwdg.size()==0) {
      for (int i=0;i<which_nucl.size();i++){
      	COGwdg.push_back(1);
      }
    } else {
      if (COGwdg.size()!=which_nucl.size()) {
	cerr<< "COGwdg must match the number of pulses, "<<which_nucl.size()<<", aborting.....\n";
	exit(1);
      }
      nu=COGwdg;
    }

    // the last one can often be avoided because pn-po is the same for all pws

    if (COGmin_N==NOTASSIGNED)
      COGmin_N=2;
    if (COGinc_N==NOTASSIGNED)
      COGinc_N=1;

    for(COGN=COGmin_N;COGN<=COGmax_N;COGN+=COGinc_N) {
      // N loop
      cout<<"================================================\n";
      cout<< "COGWHEEL N="<<COGN<<" searching...\n";
      cout.flush();
      // now choose nu_i's, optimiz. would go here
      // e.g. no need for wd. numb. if only allowed pw is zero
      // another one can be later economized by the Delta nu trick
      //
      // if ZERO included, then can blank out some pulses

      maxpossib=1;
      for (int i=0; i<which_nucl.size();i++){
	if (COGwdg[i]!=0)
	  maxpossib*=COGN;
      }

      for (int possib=0;possib<maxpossib;possib++){
	dummy=possib;
	// set nu_i
	for (int pind=0;pind<which_nucl.size();pind++){
	  if(COGwdg[pind]!=0){
	    nu[pind]=dummy%COGN;
	    dummy=dummy/COGN;
	  }
	}
	//cout <<"nu "<<nu<<endl;

	initialize_pw_idx(pw_idx,newattn);

	// loop over pathways
	paths.clear();
	do {
	  make_pw(currpw,pw_idx,newattn);
	  //cout << currpw;
	  // for cogwheel: can disregard last Delta p if starting and ending coh. ord same as
	  // for desired pw.?
	  make_diffpw(diffpw,currpw,index_to_currpw,which_nucl,phase_allequal,nucl_pw_idx);
	  // calc. phase cycling sig
	  sum=0;
	  for (int i=0;i<which_nucl.size();i++)
	    sum+=nu[i]*diffpw[i];
	  sum=sum%COGN;

	  if (sum==0){
	    paths.push_back(currpw);
	  }
	} while (increment_pw_idx(pw_idx,newattn,newattn_max)&&(paths.size()<=COGmax_pw));

	if (paths.size()<=COGmax_pw){
	  cout <<"found: COG "<<COGN<<" (";
	  for (int i=0;i<nu.size();i++){
	    cout <<nu[i];
	    if (i<nu.size()-1)
	      cout<<",";
	  }
	  int nu_r=calc_cogwheel_receiver_phase(which_nucl,nu,COGpw,COGN);
	  cout <<"; "<<nu_r<<")"<<endl;
	  cout <<"selected pathways: "<<paths.size()<<endl;
	  outputpaths(paths,COGpw);
	  cout <<"--------------------------------------------\n";
	}
      }
    };
   }
}

	  
	















