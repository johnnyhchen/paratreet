#ifndef SIMPLE_CACHEMANAGER_H_
#define SIMPLE_CACHEMANAGER_H_

#include "simple.decl.h"
#include "common.h"
#include "MultiMsg.h"
#include "Utility.h"
#include <map>
#include <unordered_map>
#include <vector>
#include "templates.h"

template <typename Data>
class CacheManager : public CBase_CacheManager<Data> {
public:
  CProxy_TreePiece<Data> tp_proxy;
#ifdef SMPCACHE
  CmiNodeLock clock, tlock;
#endif
  Node<Data>* root;
  std::map<Key, Node<Data>*> tps;
  std::map<Key, std::set<int> > curr_waiting;
  std::vector<Node<Data>*> delete_at_end;

  CacheManager() { // : root(nullptr), curr_waiting (std::map<Key, std::vector<int> >()) {}
    root = new Node<Data>(1, 0, 0, NULL, 0, 0, NULL);
    root->type = Node<Data>::Boundary;
#ifdef SMPCACHE
    clock = CmiCreateLock();
    tlock = CmiCreateLock();
#endif
  }

  void receiveTP(TPHolder<Data> tp_holderi) {
    tp_proxy = tp_holderi.tp_proxy;
  }

  ~CacheManager() {
    for (int i = 0; i < delete_at_end.size(); i++) delete delete_at_end[i];
    root->triggerFree();
    delete root;
    curr_waiting.clear();
  }
  bool connect (Node<Data>*);
  Node<Data>* findNode(Key);
  template <typename Visitor>
  void addCache(MultiMsg<Data>*);
  template <typename Visitor>
  void addCacheHelper(Particle*, int, Node<Data>*, int);
  template <typename Visitor>
  void restoreData(Key, Data);
  void addLocal(Key, Node<Data>*);
  void processNode(Node<Data>*);
  Node<Data> *getNode(Key);
  template <typename Visitor>
  void resumeTraversals(Key);
  void insertNode(Node<Data>*, bool, bool);
  void swapIn(Node<Data>*);
};

template <typename Data>
template <typename Visitor>
void CacheManager<Data>::addCache(MultiMsg<Data>* multimsg) {
  addCacheHelper<Visitor>(multimsg->particles, multimsg->n_particles, multimsg->nodes, multimsg->n_nodes);
  delete multimsg;
}

template <typename Data>
template <typename Visitor>
void CacheManager<Data>::addCacheHelper(Particle* particles, int n_particles, Node<Data>* nodes, int n_nodes) {
  Node<Data>* first_node = findNode(nodes[0].key), *first_added_node;
  if (first_node->type != Node<Data>::CachedRemote && first_node->type != Node<Data>::CachedRemoteLeaf) {
    int p_index = 0;
    for (int j = 0; j < n_nodes; j++) {
      Node<Data>* node = new Node<Data>(nodes[j]);
      if (j == 0) first_added_node = node;
      if (node->type == Node<Data>::Leaf || node->type == Node<Data>::EmptyLeaf) {
        node->type = Node<Data>::CachedRemoteLeaf;
        if (node->n_particles) node->particles = new Particle [node->n_particles];
        for (int i = 0; i < node->n_particles; i++) {
          if (p_index < n_particles) node->particles[i] = particles[p_index++];
          else CkPrintf("yikes not good\n");
        }
      }
      else if (node->type == Node<Data>::Internal) {
        node->type = Node<Data>::CachedRemote;
      }
#ifdef SMPCACHE
      if (j == 0) insertNode(node, false, false);
      else insertNode(node, false, true);
#else
      insertNode(node, false, true);
#endif
    }
  }
#ifdef SMPCACHE
  swapIn(first_added_node);
#endif
  resumeTraversals<Visitor>(nodes[0].key);
}
/*
template <typename Data>
void CacheManager<Data>::addLocal(Key key, Node<Data> *node) {
  
}*/
/*
template<typename Data>
Node<Data>* CacheManager<Data>::getNode(Key key) {
  Node<Data> *node = root;
  int depth = Utility::getDepthFromKey(key), child_idx;
  while(depth) {
    depth -= 1;
    child_idx = (key >> (depth * LOG_BRANCH_FACTOR)) % BRANCH_FACTOR;

    if(node->n_children == 0) {
        return NULL;
    }
    node = node->children[child_idx];
  }
  return node;
}*/

template <typename Data>
template <typename Visitor>
void CacheManager<Data>::resumeTraversals (Key key) {
#ifdef SMPCACHE
  CmiLock(clock);
#endif
  if (curr_waiting.count(key)) {
    std::set<int> indices = curr_waiting[key];
    for (std::set<int>::iterator it = indices.begin(); it != indices.end(); it++) {
      //CkPrintf("restoring tp %d\n", indices[j]);
      tp_proxy[*it].template goDown<Visitor>(key);
    }
    curr_waiting.erase(key);
  }
#ifdef SMPCACHE
  CmiUnlock(clock);
#endif
}

template <typename Data>
template <typename Visitor>
void CacheManager<Data>::restoreData(Key key, Data di) {
  if (findNode(key)->type != Node<Data>::CachedBoundary) {
    //CkPrintf("restoring %d's data on pe %d\n", key, CkMyPe());
    Node<Data>* node = new Node<Data>();
    node->data = di;
    node->key = key;
    node->type = Node<Data>::CachedBoundary;
    node->n_children = 8; 
    insertNode(node, true, true);
  }
  resumeTraversals<Visitor>(key);
}

template <typename Data>
void CacheManager<Data>::swapIn(Node<Data>* to_swap) {
  Node<Data>* copy = findNode(to_swap->key);
  Node<Data>* parent = copy->parent;
  delete_at_end.push_back(copy);
  parent->children[to_swap->key % 8] = to_swap;
  to_swap->parent = copy->parent;
}

template <typename Data>
bool CacheManager<Data>::connect(Node<Data>* node) {
  Node<Data>* copy = findNode(node->key);
  if (copy == NULL) {
#ifdef SMPCACHE
    CmiLock(tlock);
#endif
    tps.insert(std::make_pair(node->key, node));
#ifdef SMPCACHE
    CmiUnlock(tlock);
#endif
    return false;
  }
  CkPrintf("didnt add to tps\n");
  // need to lock this? no i dont i dont think
  swapIn(node);
  return true;
}

template <typename Data>
void CacheManager<Data>::insertNode(Node<Data>* node, bool above_tp, bool should_swap) {
  //CkPrintf("start %d\n", CkMyPe());
  for (int i = 0; i < node->n_children; i++) {
    Node<Data>* new_child;
#ifdef SMPCACHE
    CmiLock(tlock);
#endif
    if (tps.count(node->key * 8 + i)) {
      new_child = tps[node->key * 8 + i];
      tps.erase(new_child->key);
    } 
    else {
      new_child = new Node<Data> (node->key * 8 + i, node->depth+1, 0, NULL, 0, 0, node);
      new_child->type = (above_tp) ? Node<Data>::RemoteAboveTPKey : Node<Data>::Remote;
      new_child->tp_index = node->tp_index;
    }
#ifdef SMPCACHE
    CmiUnlock(tlock);
#endif
    node->children.push_back(new_child);
  }
  if (should_swap) {
    if (node->key == 1) {
      delete_at_end.push_back(root);
      root = node;
    }
    else swapIn(node);
  }
}

/*  CmiLock(rlock);
  Node<Data> *curr = root, *parent = NULL;
  Key node_key = node->key, child_key;
  int depth = Utility::getDepthFromKey(node_key), child_idx;

  while(depth) {
    depth -= 1;
    child_idx = (node_key >> (depth * LOG_BRANCH_FACTOR)) % BRANCH_FACTOR;

    if(curr->n_children == 0) {
      curr->children.resize(BRANCH_FACTOR, NULL);
      curr->n_children = BRANCH_FACTOR;
    }

    if(curr->children[child_idx] == NULL) {
      child_key = (curr->key) << LOG_BRANCH_FACTOR + child_idx;
      if(child_key == node_key) {
        curr->children[child_idx] = node;
        node->parent = curr;
        return;
      }
      curr->children[child_idx] = new Node<Data>(child_key, curr->depth + 1, 0, NULL, 0, 0, curr);
    }
    parent = curr;
    curr = curr->children[child_idx];
  }

  if(curr->type != 0) {
     //non Invalid node - possible?
     return;
  }

  if(parent) {
    parent->children[child_idx] = node;
  } else {
    root = node;
  }
  node->parent = parent;
  node->depth = curr->depth;
  if(curr->n_children) {
    node->children.assign(curr->children.begin(), curr->children.end());
    curr->children.clear();
  }
  //anything else we need from Invalid node?
  //particles?
  delete curr; //it shouldn't have any particles? Invalid node
  CmiUnlock(rlock);
} */

template <typename Data>
Node<Data>* CacheManager<Data>::findNode(Key key) {
  std::vector<int> remainders; // sorry lol
  Key temp = key;
  while (temp >= 8) {
    remainders.push_back(temp % 8);
    temp /= 8;
  }
  Node<Data>* node = root;
  for (int i = remainders.size()-1; i >= 0; i--) {
    if (remainders[i] < node->children.size()) node = node->children[remainders[i]];
    else return NULL;
  }
  return node;
}

#endif //SIMPLE_CACHEMANAGER_H_
