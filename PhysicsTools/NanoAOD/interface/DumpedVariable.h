#ifndef PhysicsTools_NanoAOD_DumpedVariable_h
#define PhysicsTools_NanoAOD_DumpedVariable_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <string>

// Base class for dumped variables
class VariableBase {
public:
  VariableBase(const std::string &aname, const edm::ParameterSet &cfg)
      : name_(aname),
        doc_(cfg.getParameter<std::string>("doc")),
        precision_(cfg.existsAs<int>("precision") ? cfg.getParameter<int>("precision")
                                                  : (cfg.existsAs<std::string>("precision") ? -2 : -1)) {}
  virtual ~VariableBase() {}
  const std::string &name() const { return name_; }

protected:
  std::string name_, doc_;
  int precision_;
};

// Object member variables and methods
template <typename ObjType>
class Variable : public VariableBase {
public:
  Variable(const std::string &aname, const edm::ParameterSet &cfg) : VariableBase(aname, cfg) {}
  virtual void fill(std::vector<const ObjType *> &selobjs, nanoaod::FlatTable &out) const = 0;
};

template <typename ObjType, typename StringFunctor, typename ValType>
class FuncVariable : public Variable<ObjType> {
public:
  FuncVariable(const std::string &aname, const edm::ParameterSet &cfg)
      : Variable<ObjType>(aname, cfg),
        func_(cfg.getParameter<std::string>("expr"), true),
        precisionFunc_(cfg.existsAs<std::string>("precision") ? cfg.getParameter<std::string>("precision") : "23",
                       true) {}
  ~FuncVariable() override {}
  void fill(std::vector<const ObjType *> &selobjs, nanoaod::FlatTable &out) const override {
    std::vector<ValType> vals(selobjs.size());
    for (unsigned int i = 0, n = vals.size(); i < n; ++i) {
      ValType val = func_(*selobjs[i]);
      if constexpr (std::is_same<ValType, float>()) {
        if (this->precision_ == -2) {
          auto prec = precisionFunc_(*selobjs[i]);
          vals[i] = prec > 0 ? MiniFloatConverter::reduceMantissaToNbitsRounding(val, prec) : val;
        } else
          vals[i] = val;
      } else {
        vals[i] = val;
      }
    }
    out.template addColumn<ValType>(this->name_, vals, this->doc_, this->precision_);
  }

protected:
  StringFunctor func_;
  StringFunctor precisionFunc_;
};

// External variables: i.e. variables that are not member or methods of the object
template <typename ObjType>
class ExtVariable : public VariableBase {
public:
  ExtVariable(const std::string &aname, const edm::ParameterSet &cfg) : VariableBase(aname, cfg) {}
  virtual void fill(const edm::Event &iEvent,
                    std::vector<edm::Ptr<ObjType>> selptrs,
                    nanoaod::FlatTable &out) const = 0;
};
template <typename ObjType, typename TIn, typename ValType = TIn>
class ValueMapVariable : public ExtVariable<ObjType> {
public:
  ValueMapVariable(const std::string &aname,
                   const edm::ParameterSet &cfg,
                   edm::ConsumesCollector &&cc,
                   bool skipNonExistingSrc = false)
      : ExtVariable<ObjType>(aname, cfg),
        skipNonExistingSrc_(skipNonExistingSrc),
        token_(cc.consumes<edm::ValueMap<TIn>>(cfg.getParameter<edm::InputTag>("src"))) {}
  void fill(const edm::Event &iEvent, std::vector<edm::Ptr<ObjType>> selptrs, nanoaod::FlatTable &out) const override {
    edm::Handle<edm::ValueMap<TIn>> vmap;
    iEvent.getByToken(token_, vmap);
    std::vector<ValType> vals;
    if (vmap.isValid() || !skipNonExistingSrc_) {
      vals.resize(selptrs.size());
      for (unsigned int i = 0, n = vals.size(); i < n; ++i) {
        vals[i] = (*vmap)[selptrs[i]];
      }
    }
    out.template addColumn<ValType>(this->name_, vals, this->doc_, this->precision_);
  }

protected:
  const bool skipNonExistingSrc_;
  edm::EDGetTokenT<edm::ValueMap<TIn>> token_;
};

#endif
