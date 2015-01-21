//New class
#ifndef TREEVARIABLES
#define TREEVARIABLES
/////
//   Headers
/////
#include <TTree.h>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
/////
//   Constants
/////
#define DEF_SIZE1D 25
#define DEF_VAL_INT -9999
#define DEF_VAL_FLOAT -9999.0f
#define DEF_VAL_DOUBLE -9999.0d
#define FLOAT_EPS 0.0000001f
#define DOUBLE_EPS 0.0000001d
/////
//   Functions
/////
#define INIT_1DARRAY(x,n,y) for(int i=0;i<n;i++) {x[i]=y;}
#define INIT_2DARRAY(x,n,m,y) for(int i=0;i<n;i++) { for(int j=0;j<m;j++) { x[i][j]=y; } }
inline bool is_undef(int x) { return x==DEF_VAL_INT; };
inline bool is_undef(float x) { return fabs(x-DEF_VAL_FLOAT) < FLOAT_EPS; };
inline bool is_undef(double x) { return fabs(x-DEF_VAL_DOUBLE) < DOUBLE_EPS; }
/////
//   Class declaration
/////
class CTree{
public:
 CTree(TTree* _tree) { tree = _tree; };
 TTree* tree;
 /////
 //   Helper functions for accessing branches
 /////
 template <typename T>
 T get_address(const std::string name){
  auto* br = tree->GetBranch(name.c_str());
  if(br==0){
   std::cerr << "ERROR: get_address CTree " << "branch " << name << " does not exist" << std::endl;
   throw std::exception();
  }
  auto* p = br->GetAddress();
  return reinterpret_cast<T>(p);
 }
 /////
 //   Declare variables
 /////
 //Event
 int evt_id, evt_json, evt_lumi, evt_run;
 //Test vars
 double prova;
 int    provarrsize;
 double provarr[DEF_SIZE1D];
 /////
 //   Initialise
 /////
 void loop_initialize(void){
  evt_id   = DEF_VAL_INT;
  evt_json = DEF_VAL_INT;
  evt_lumi = DEF_VAL_INT;
  evt_run  = DEF_VAL_INT;
  prova       = DEF_VAL_DOUBLE;
  provarrsize = DEF_SIZE1D;
  INIT_1DARRAY(provarr,DEF_SIZE1D,DEF_VAL_DOUBLE); 
 } 
 /////
 //   Set branches
 /////
 void make_branches(void){
  tree->Branch("evt_id", &evt_id, "evt_id/I");
  tree->Branch("evt_json", &evt_json, "evt_json/I");
  tree->Branch("evt_lumi", &evt_lumi, "evt_lumi/I");
  tree->Branch("evt_run", &evt_run, "evt_run/I");
  tree->Branch("prova", &prova, "prova/D");
  tree->Branch("provarrsize", &provarrsize, "provarrsize/I");
  tree->Branch("provarr", provarr, "provarr[provarrsize]/D");
 }
 /////
 //   Set branch address
 /////
 //Connects the branches of an existing TTree to variables used when loading the file
 void set_branch_addresses(void){
  tree->SetBranchAddress("evt_id", &evt_id);
  tree->SetBranchAddress("evt_json", &evt_json);
  tree->SetBranchAddress("evt_lumi", &evt_lumi);
  tree->SetBranchAddress("evt_run", &evt_run);
  tree->SetBranchAddress("prova", &prova);
  tree->SetBranchAddress("provarrsize", &provarrsize);
  tree->SetBranchAddress("provarr", &provarr);
 }
};
#endif
