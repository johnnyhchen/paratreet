#ifndef SIMPLE_TREEELEMENT_H_
#define SIMPLE_TREEELEMENT_H_

#include "simple.decl.h"

template<typename Visitor, typename Data>
class TreeElement : public CBase_TreeElement<Visitor, Data> {
private:
  Visitor v;
  Data d;
  int wait_count;
public:
  TreeElement();
  void receiveData (Data, bool);
};

extern CProxy_Main mainProxy;

template<typename Visitor, typename Data>
TreeElement<Visitor, Data>::TreeElement() :
  v (Visitor()), d (Data()), wait_count (-1) {
  //CkPrintf("%d\n", thisIndex);
}
template<typename Visitor, typename Data>
void TreeElement<Visitor, Data>::receiveData (Data di, bool if_leafi) {
  if (wait_count == -1) wait_count = (if_leafi) ? 1 : 8;
  d = d + di;
  wait_count--;
  if (wait_count == 0) {
    v.node(mainProxy, this->thisProxy, d, this->thisIndex);
  }
}

#define CK_TEMPLATES_ONLY
#include "simple.def.h"
#undef CK_TEMPLATES_ONLY

#endif // SIMPLE_TREEELEMENT_H_
