#include "myVtxObj.h"
#ifdef __VTX_HADRON_COLLECTION__


using std::vector;
using std::string;

myVtx::hadronAtt::hadronAtt (const int id, 
    const int ipt, const float eta, const float phi):
    id_(id), ipt_(ipt), eta_(eta), phi_(phi) 
{ }

myVtx::hadronAtt::hadronAtt (const myVtx::hadronAtt& h)
{
  id_ = h.id();
  ipt_ = h.ipt();
  eta_ = h.eta();
  phi_ = h.phi();
}

myVtx::hadronAtt::~hadronAtt() { }

void myVtx::hadronAtt::setAtt(const int id, 
    const int ipt, const float eta, const float phi)
{
  id_ = id;
  ipt_ = ipt;
  eta_ = eta;
  phi_ = phi;
}

void myVtx::hadronAtt::clear()
{
  id_ = -99;
  ipt_ = -99;
  eta_ = -99;
  phi_ = -99;
}

myVtx::vtxAtt::vtxAtt(const int id, 
    const int iz, const int itrk):
  id_(id), iz_(iz), itrk_(itrk)
{ }

myVtx::vtxAtt::~vtxAtt() { }

void myVtx::vtxAtt::addVecHadron(const std::vector<myVtx::hadronAtt>& vec_h, const string& type)
{
  if(type == "trig") {
    vec_trig_.insert( vec_trig_.end(), vec_h.begin(), vec_h.end() );
  }
  else if(type == "assc"){
    vec_assc_.insert( vec_assc_.end(), vec_h.begin(), vec_h.end() );
  }
  return;
}

#endif
