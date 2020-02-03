#ifndef __VTX_HADRON_COLLECTION__
#define __VTX_HADRON_COLLECTION__

#include <cstddef>
#include <vector>
#include <string>

namespace myVtx{

  class hadronAtt {
    public: 
      explicit hadronAtt(const int, const int, 
        const float, const float);
      hadronAtt(const hadronAtt&);
      void setAtt(const int, const int, 
        const float, const float);
      void clear();
      ~hadronAtt();

      int id() const { return id_; }
      int ipt() const { return ipt_; }
      float eta() const { return eta_; }
      float phi() const { return phi_; }

    private:
      // I do not know how to correctly define operator=
      int    id_;
      int    ipt_;
      float eta_;
      float phi_;
  };

  class vtxAtt {
    public:
      explicit vtxAtt(const int, const int, const int);
      ~vtxAtt();
      const std::vector<myVtx::hadronAtt>& vec_trig() { return vec_trig_; }
      const std::vector<myVtx::hadronAtt>& vec_assc() { return vec_assc_; }
      int id() const { return id_; }
      int iz() const { return iz_; }
      int itrkBin() const { return itrk_; }

      void addVecHadron(const std::vector<hadronAtt>&, const std::string&);
    private:
      int   id_;
      int   iz_;
      int   itrk_;
      std::vector<hadronAtt> vec_trig_;
      std::vector<hadronAtt> vec_assc_;
  };

  template <typename T>
    int getIdx(const std::vector<T>& v, const T& val)
    {
      int size = v.size() - 1;
      for(int iv=0; iv<size; iv++){
        if(val < v[iv+1] && val > v[iv]) return iv;
      }
    return -1;
    }
};

#endif
