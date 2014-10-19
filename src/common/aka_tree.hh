/**
 * @file   aka_tree.hh
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Tue May 13 2014
 *
 * @brief  bounding volume hierarchies
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_TREE_HH__
#define __AKANTU_TREE_HH__

#include <cassert>
#include <queue>
#include <set>
#include <iomanip>
//#include <limits>

#include "aka_ball.hh"
#include "aka_bounding_box.hh"
#include "aka_array.hh"

//#define DEBUG_TREE 1


__BEGIN_AKANTU__


using std::cout;
using std::endl;


struct Tree_node_base
{
  typedef Tree_node_base* base_ptr;
  typedef const Tree_node_base* const_base_ptr;
  
  base_ptr parent_;
  base_ptr left_;
  base_ptr right_;
  
  static base_ptr _minimum(base_ptr x) {
    while (x->left_ != 0) x = x->left_;
    return x;
  }
  
  static const_base_ptr _minimum(const_base_ptr x) {
    while (x->left_ != 0) x = x->left_;
    return x;
  }
  
  static base_ptr _maximum(base_ptr x) {
    while (x->right_ != 0) x = x->right_;
    return x;
  }
  
  static const_base_ptr _maximum(const_base_ptr x) {
    while (x->right_ != 0) x = x->right_;
    return x;
  }
};

template<typename T>
struct Tree_node : public Tree_node_base
{
  typedef T value_tpe;
  typedef Tree_node<value_tpe>* link_type;
  
  value_tpe value_;
  
  template<typename... Args>
  Tree_node(Args&&... args)
	: Tree_node_base(), value_(std::forward<Args>(args)...) { }
};


// iterators

static Tree_node_base* local_leaf_increment(Tree_node_base* x) throw () {
  
  Tree_node_base* y = x->parent_;
  if (y->parent_ == x)
  return y;
  while (x == y->right_) {
    x = y;
    y = y->parent_;
  }
  if (y->parent_ == x)
  return y;
  if (x->right_ != y)
  x = y;
  if (x->right_ != 0)  {
    x = x->right_;
    while (x->left_ != 0)
    x = x->left_;
  }
  return x;
}

template <class tree_pointer>
tree_pointer* _tree_leaf_increment(tree_pointer* x) throw ()
{ return local_leaf_increment(x); }

template <class tree_pointer>
const tree_pointer* _tree_leaf_increment(const tree_pointer* x) throw ()
{ return local_leaf_increment(const_cast<tree_pointer*>(x)); }


template<typename _Tp>
struct Tree_leaf_iterator
{
  typedef _Tp  value_type;
  typedef _Tp& reference;
  typedef _Tp* pointer;
  
  typedef std::bidirectional_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;
  
  typedef Tree_leaf_iterator<_Tp> self_type;
  typedef Tree_node_base::base_ptr base_ptr;
  typedef Tree_node<_Tp>* link_type;
  
  Tree_leaf_iterator() : node_() { }
  
  explicit Tree_leaf_iterator(link_type x) : node_(x) { }
  
  reference operator*() const
  { return static_cast<link_type>(node_)->value_; }
  
  pointer operator->() const
  { return std::__addressof(static_cast<link_type>(node_)->value_); }
  
  self_type& operator++() {
    node_ = _tree_leaf_increment(node_);
    return *this;
  }
  
  self_type operator++(int) {
    self_type tmp = *this;
    node_ = _tree_leaf_increment(node_);
    return tmp;
  }
  
  bool operator==(const self_type& x) const
  { return node_ == x.node_; }
  
  bool operator!=(const self_type& x) const
  { return node_ != x.node_; }
  
  bool operator<(const self_type& x) const
  { return node_ < x.node_; }
  
  base_ptr node_;
};


template<typename _Tp>
struct Tree_leaf_const_iterator
{
  typedef _Tp        value_type;
  typedef const _Tp& reference;
  typedef const _Tp* pointer;
  
  typedef Tree_leaf_iterator<_Tp> iterator;
  
  typedef std::bidirectional_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;
  
  typedef Tree_leaf_const_iterator<_Tp> self_type;
  typedef Tree_node_base::const_base_ptr base_ptr;
  typedef const Tree_node<_Tp>* link_type;
  
  Tree_leaf_const_iterator() : node_() { }
  
  explicit Tree_leaf_const_iterator(link_type x) : node_(x) { }
  
  Tree_leaf_const_iterator(const iterator& it) : node_(it.node_) { }
  
  iterator _const_cast() const
  { return iterator(static_cast<typename iterator::link_type>(const_cast<typename iterator::base_ptr>(node_))); }
  
  reference operator*() const
  { return static_cast<link_type>(node_)->value_; }
  
  pointer operator->() const
  { return std::__addressof(static_cast<link_type>
                            (node_)->value_); }
  
  self_type& operator++() {
    node_ = _tree_leaf_increment(node_);
    return *this;
  }
  
  self_type operator++(int) {
    self_type tmp = *this;
    node_ = _tree_leaf_increment(node_);
    return tmp;
  }
  
  bool operator==(const self_type& x) const
  { return node_ == x.node_; }
  
  bool operator!=(const self_type& x) const
  { return node_ != x.node_; }
  
  bool operator<(const self_type& x) const
  { return node_ < x.node_; }
  
  base_ptr node_;
};


template<typename _Val>
inline bool operator==(const Tree_leaf_iterator<_Val>& x,
                       const Tree_leaf_const_iterator<_Val>& y)
{ return x.node_ == y.node_; }

template<typename _Val>
inline bool operator!=(const Tree_leaf_iterator<_Val>& x,
                       const Tree_leaf_const_iterator<_Val>& y)
{ return x.node_ != y.node_; }



template <class tree_pointer>
static tree_pointer*
local_tree_increment(tree_pointer* x) throw ()
{
  if (x->right_ != 0)
  {
    x = x->right_;
    while (x->left_ != 0)
    x = x->left_;
  }
  else
  {
    tree_pointer* y = x->parent_;
    while (x == y->right_)
    {
      x = y;
      y = y->parent_;
    }
    if (x->right_ != y)
    x = y;
  }
  return x;
}


template <class tree_pointer>
tree_pointer* _tree_increment(tree_pointer* x) throw ()
{ return local_tree_increment(x); }

template <class tree_pointer>
const tree_pointer*
_tree_increment(const tree_pointer* x) throw ()
{ return local_tree_increment(const_cast<tree_pointer*>(x)); }

template <class tree_pointer>
static tree_pointer*
local_tree_decrement(tree_pointer* x) throw () {
  if (x->parent_->parent_ == x)
  x = x->right_;
  else if (x->left_ != 0) {
    Tree_node_base* y = x->left_;
    while (y->right_ != 0)
    y = y->right_;
    x = y;
  }
  else {
    tree_pointer* y = x->parent_;
    while (x == y->left_) {
      x = y;
      y = y->parent_;
    }
    x = y;
  }
  return x;
}

template <class tree_pointer>
tree_pointer* _tree_decrement(tree_pointer* x) throw ()
{ return local_tree_decrement(x); }

template <class tree_pointer>
const tree_pointer* _tree_decrement(const tree_pointer* x) throw ()
{ return local_tree_decrement(const_cast<tree_pointer*>(x)); }


// iterators

template<typename _Tp>
struct Tree_iterator
{
  typedef _Tp  value_type;
  typedef _Tp& reference;
  typedef _Tp* pointer;
  
  typedef std::bidirectional_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;
  
  typedef Tree_iterator<_Tp> self_type;
  typedef Tree_node_base::base_ptr base_ptr;
  typedef Tree_node<_Tp>* link_type;
  
  typedef Tree_leaf_iterator<value_type> leaf_iterator;
  
  
  Tree_iterator() : node_() { }
  
  explicit Tree_iterator(link_type x) : node_(x) { }
  
  Tree_iterator(const leaf_iterator& it) : node_(it.node_) { }
  
  reference operator*() const
  { return static_cast<link_type>(node_)->value_; }
  
  pointer operator->() const
  { return std::__addressof(static_cast<link_type>(node_)->value_); }
  
  self_type& operator++() {
    node_ = _tree_increment(node_);
    return *this;
  }
  
  self_type operator++(int) {
    self_type tmp = *this;
    node_ = _tree_increment(node_);
    return tmp;
  }
  
  self_type& operator--() {
    node_ = _tree_decrement(node_);
    return *this;
  }
  
  self_type operator--(int) {
    self_type tmp = *this;
    node_ = _tree_decrement(node_);
    return tmp;
  }
  
  
  self_type left()
  { return self_type(static_cast<link_type>(node_->left_)); }
  
  self_type right()
  { return self_type(static_cast<link_type>(node_->right_)); }
  
  self_type parent()
  { return self_type(static_cast<link_type>(node_->parent_)); }
  
  bool is_leaf() {
    assert(node_ != nullptr);
    return node_->left_ == nullptr && node_->right_ == nullptr;
  }
  
  
  bool operator==(const self_type& x) const
  { return node_ == x.node_; }
  
  bool operator!=(const self_type& x) const
  { return node_ != x.node_; }
  
  bool operator<(const self_type& x) const
  { return node_ < x.node_; }
  
  base_ptr node_;
};

template<typename _Tp>
struct Tree_const_iterator
{
  typedef _Tp        value_type;
  typedef const _Tp& reference;
  typedef const _Tp* pointer;
  
  typedef Tree_iterator<_Tp> iterator;
  
  typedef std::bidirectional_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;
  
  typedef Tree_const_iterator<_Tp> self_type;
  typedef Tree_node_base::const_base_ptr base_ptr;
  typedef const Tree_node<_Tp>* link_type;
  
  Tree_const_iterator() : node_() { }
  
  explicit Tree_const_iterator(link_type x) : node_(x) { }
  
  Tree_const_iterator(const iterator& it) : node_(it.node_) { }
  
  iterator _const_cast() const
  { return iterator(static_cast<typename iterator::link_type>(const_cast<typename iterator::base_ptr>(node_))); }
  
  reference operator*() const
  { return static_cast<link_type>(node_)->value_; }
  
  pointer operator->() const
  { return std::__addressof(static_cast<link_type>
                            (node_)->value_); }
  
  self_type& operator++() {
    node_ = _tree_increment(node_);
    return *this;
  }
  
  self_type operator++(int) {
    self_type tmp = *this;
    node_ = _tree_increment(node_);
    return tmp;
  }
  
  self_type& operator--() {
    node_ = _tree_decrement(node_);
    return *this;
  }
  
  self_type operator--(int) {
    self_type tmp = *this;
    node_ = _tree_decrement(node_);
    return tmp;
  }
  
  self_type left()
  { return self_type(static_cast<link_type>(node_->left_)); }
  
  self_type right()
  { return self_type(static_cast<link_type>(node_->right_)); }
  
  self_type parent()
  { return self_type(static_cast<link_type>(node_->parent_)); }
  
  bool is_leaf() {
    assert(node_ != nullptr);
    return node_->left_ == nullptr && node_->right_ == nullptr;
  }
  
  
  bool operator==(const self_type& x) const
  { return node_ == x.node_; }
  
  bool operator!=(const self_type& x) const
  { return node_ != x.node_; }
  
  bool operator<(const self_type& x) const
  { return node_ < x.node_; }
  
  base_ptr node_;
};

template<typename _Val>
inline bool operator==(const Tree_iterator<_Val>& x,
                       const Tree_const_iterator<_Val>& y)
{ return x.node_ == y.node_; }

template<typename _Val>
inline bool operator!=(const Tree_iterator<_Val>& x,
                       const Tree_const_iterator<_Val>& y)
{ return x.node_ != y.node_; }




// traversal
template <class iterator_type, class functor_type>
void inorder(iterator_type it, functor_type& fn) {
  if (it.node_ == nullptr)
  return;
  inorder(it.left(), fn);
  fn(it);
  inorder(it.right(), fn);
}

template <class iterator_type, class functor_type>
void preorder(iterator_type it, functor_type& fn) {
  if (it.node_ == nullptr)
  return;
  fn(it);
  preorder(it.left(), fn);
  preorder(it.right(), fn);
}

template <class iterator_type, class functor_type>
void postorder(iterator_type it, functor_type& fn) {
  if (it.node_ == nullptr)
  return;
  postorder(it.left(), fn);
  postorder(it.right(), fn);
  fn(it);
}



//! Tree class
template<typename T, template <class> class _Cost, typename _Alloc = std::allocator<T> >
class Tree
{
  typedef typename _Alloc::template rebind<Tree_node<T> >::other
  _Node_allocator;
  
  protected:
  typedef Tree_node_base* base_ptr;
  typedef const Tree_node_base* const_base_ptr;
  
  public:
  
  typedef T value_type;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef Tree_node<T>* link_type;
  typedef const Tree_node<T>* const_link_type;
  typedef size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef _Cost<value_type> cost_functor;
  typedef _Alloc allocator_type;
  
  _Node_allocator& _get_Node_allocator()
  { return *static_cast<_Node_allocator*>(&this->impl_); }
  
  const _Node_allocator& _get_Node_allocator() const
  { return *static_cast<const _Node_allocator*>(&this->impl_); }
  
  allocator_type
  get_allocator() const
  { return allocator_type(_get_Node_allocator()); }
  
  private:
  
  link_type _get_node()
  { return impl_._Node_allocator::allocate(1); }
  
  void _put_node(link_type p)
  { impl_._Node_allocator::deallocate(p, 1); }
  
  void _destroy_node(link_type p) {
    _get_Node_allocator().destroy(p);
    _put_node(p);
  }
  
  link_type _clone_node(const_link_type x) {
    link_type tmp = _create_node(x->value_);
    tmp->left_ = 0;
    tmp->right_ = 0;
    return tmp;
  }
  
  template<typename... _Args>
  link_type _create_node(_Args&&... args) {
    
    link_type tmp = _get_node();
    try {
      _get_Node_allocator().construct(tmp, std::forward<_Args>(args)...);
      tmp->left_ = 0;
      tmp->right_ = 0;
      tmp->parent_ = 0;
    }
    catch(...) {
      _put_node(tmp);
      throw;
    }
    return tmp;
	}
  
  private:
  
  template<typename _Cost_functor>
  struct Tree_impl : public _Node_allocator
  {
    _Cost_functor cost_;
    Tree_node_base header_;
    size_type node_count_; // Keeps track of size of tree.
    
    Tree_impl()
    : _Node_allocator(), cost_(), header_(),
    node_count_(0)
    { _M_initialize(); }
    
    Tree_impl(const _Cost_functor& comp, const _Node_allocator& a)
    : _Node_allocator(a), cost_(comp), header_(),
    node_count_(0)
    { _M_initialize(); }
    
    Tree_impl(const _Cost_functor& comp, _Node_allocator&& a)
    : _Node_allocator(std::move(a)), cost_(comp),
    header_(), node_count_(0)
    { _M_initialize(); }
    
    private:
    void
    _M_initialize()
    {
      this->header_.parent_ = 0;
      this->header_.left_ = &this->header_;
      this->header_.right_ = &this->header_;
    }
	};
  
  private:
  
  Tree_impl<cost_functor> impl_;          //!< Tree implementation instance
  
  base_ptr& _root()
  { return this->impl_.header_.parent_; }
  
  const_base_ptr _root() const
  { return this->impl_.header_.parent_; }
  
  base_ptr& _leftmost()
  { return this->impl_.header_.left_; }
  
  const_base_ptr _leftmost() const
  { return this->impl_.header_.left_; }
  
  base_ptr& _rightmost()
  { return this->impl_.header_.right_; }
  
  const_base_ptr _rightmost() const
  { return this->impl_.header_.right_; }
  
  link_type _begin()
  { return static_cast<link_type>(this->impl_.header_.parent_); }
  
  const_link_type _begin() const
  { return static_cast<const_link_type>(this->impl_.header_.parent_); }
  
  link_type _end()
  { return static_cast<link_type>(&this->impl_.header_); }
  
  const_link_type _end() const
  { return static_cast<const_link_type>(&this->impl_.header_); }
  
  static reference _value(link_type x)
  { return x->value_; }
  
  static const_reference _value(const_link_type x)
  { return x->value_; }
  
  static reference _value(base_ptr x)
  { return static_cast<link_type>(x)->value_; }
  
  static const_reference _value(const_base_ptr x)
  { return static_cast<const_link_type>(x)->value_; }
  
  static link_type _left(base_ptr x)
  { return static_cast<link_type>(x->left_); }
  
  static const_link_type _left(const_base_ptr x)
  { return static_cast<const_link_type>(x->left_); }
  
  static link_type _right(base_ptr x)
  { return static_cast<link_type>(x->right_); }
  
  static const_link_type _right(const_base_ptr x)
  { return static_cast<const_link_type>(x->right_); }
  
  static link_type _parent(base_ptr x)
  { return static_cast<link_type>(x->parent_); }
  
  static const_link_type _parent(const_base_ptr x)
  { return static_cast<const_link_type>(x->parent_); }
  
  
  void _erase(link_type x) {
    while (x != 0)
    {
#ifdef DEBUG_TREE
      cout<<"pointer x -> "<<x<<endl;
      cout<<"calling erase on right -> "<<_right(x)<<endl;
      cout<<"destroying node -> "<<_value(x)<<endl;
      cout<<"setting x to left -> "<<_left(x)<<endl;
#endif
      _erase(_right(x));
      link_type y = _left(x);
      _destroy_node(x);
      --impl_.node_count_; // decrement nodes
      x = y;
    }
  }
  
  link_type _copy(const_link_type x, link_type p) {
    
    link_type top = _clone_node(x);
    top->parent_ = p;
    
    try
    {
      if (x->right_)
      top->right_ = _copy(_right(x), top);
      p = top;
      x = _left(x);
      
      while (x != 0)
      {
        link_type y = _clone_node(x);
        p->left_ = y;
        y->parent_ = p;
        if (x->right_)
        y->right_ = _copy(_right(x), y);
        p = y;
        x = _left(x);
      }
    }
    catch(...)
    {
      _erase(top);
      throw;
    }
    return top;
  }
  
  void _check_integrity(link_type x) {
    
    while (x != 0) {
      
      link_type r = _right(x);
      if (r)
      if (r->parent_ != x) {
        cout<<"*** ERROR *** Integrity check failed"<<endl;
        cout<<"Right child parent mismatch"<<endl;
        exit(1);
      }
      link_type l = _left(x);
      if (l)
      if (l->parent_ != x) {
        cout<<"*** ERROR *** Integrity check failed"<<endl;
        cout<<"Left child parent mismatch"<<endl;
        exit(1);
      }
      _check_integrity(r);
      x = _left(x);
    }
  }
  
  
  public:
  
  typedef Tree_iterator<value_type> iterator;
  typedef Tree_const_iterator<value_type> const_iterator;
  
  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  
  typedef Tree_leaf_iterator<value_type> leaf_iterator;
  typedef Tree_leaf_const_iterator<value_type> const_leaf_iterator;
  
  
  
  // allocation/deallocation
  Tree() {}
  
  Tree(const cost_functor& cost, const allocator_type& a = allocator_type())
  : impl_(cost, _Node_allocator(a)) { }
  
  Tree(const Tree& x) : impl_(x.impl_.cost_, x._get_Node_allocator())
  {
    if (x._root() != 0)
    {
      _root() = _copy(x._begin(), _end());
      _leftmost() = Tree_node_base::_minimum(_root());
      _rightmost() = Tree_node_base::_maximum(_root());
      impl_.node_count_ = x.impl_.node_count_;
    }
  }
  
  //! construct a tree from link type, defined by the user in any way
  Tree(link_type x) {
    
    // integrity check
    _check_integrity(x);
    if (x) _insert(0,x);
  }
  
  
  ~Tree() noexcept
  { _erase(_begin()); }
  
  const cost_functor& cost() const
  { return impl_.cost_; }
  
  void clear() noexcept
  {
    _erase(_begin());
    _leftmost() = _end();
    _root() = 0;
    _rightmost() = _end();
    // assert count == 0, _erase should take care of decrementing the node count
    assert(impl_.node_count_ == 0);
    //        impl_.node_count_ = 0;
  }
  
  iterator begin()
  { return iterator(static_cast<link_type>(this->impl_.header_.left_)); }
  
  const_iterator begin() const
  { return const_iterator(static_cast<const_link_type>(this->impl_.header_.left_)); }
  
  iterator end()
  { return iterator(static_cast<link_type>(&this->impl_.header_)); }
  
  const_iterator end() const
  { return const_iterator(static_cast<const_link_type>(&this->impl_.header_)); }
  
  
  reverse_iterator rbegin()
  { return reverse_iterator(end()); }
  
  const_reverse_iterator rbegin() const
  { return const_reverse_iterator(end()); }
  
  reverse_iterator rend()
  { return reverse_iterator(begin()); }
  
  const_reverse_iterator rend() const
  { return const_reverse_iterator(begin()); }
  
  
  leaf_iterator leaves_begin()
  { return leaf_iterator(static_cast<link_type>(impl_.header_.left_)); }
  
  const_leaf_iterator leaves_begin() const
  { return const_leaf_iterator(static_cast<const_link_type>(impl_.header_.left_)); }
  
  leaf_iterator leaves_end()
  { return leaf_iterator(static_cast<link_type>(&impl_.header_)); }
  
  const_leaf_iterator leaves_end() const
  { return const_leaf_iterator(static_cast<const_link_type>(&impl_.header_)); }
  
  iterator root()
  { return iterator(static_cast<link_type>(impl_.header_.parent_)); }
  
  const_iterator root() const
  { return const_iterator(static_cast<const_link_type>(impl_.header_.parent_)); }
  
  bool empty() const
  { return impl_.node_count_ == 0; }
  
  size_type size() const
  { return impl_.node_count_; }
  
  size_type max_size() const
  { return _get_Node_allocator().max_size(); }
  
  
  std::tuple<size_t, size_t, size_t> height() const {
    
    size_t count = 0;
    size_t min = std::numeric_limits<size_t>::max();
    size_t max = 0;
    size_t h = 0;
    for (const_leaf_iterator it = leaves_begin(); it != leaves_end(); ++it) {
      const_base_ptr x = it.node_;
      size_t h1 = 0;
      while (x->parent_ != _end()) {
        ++h1;
        x = x->parent_;
      }
      min = std::min(min, h1);
      max = std::max(max, h1);
      ++count;
      h += h1;
    }
    
    return std::make_tuple(h/count, min, max);
  }
  
  
  private:
  
  void _add_level(std::queue<link_type>& q, base_ptr x, size_t h, size_t level) {
    
    if (h < level) {
      if (x->left_ != 0) _add_level(q, x->left_, h+1, level);
      if (x->right_ != 0) _add_level(q, x->right_, h+1, level);
    }
    else if (h == level)
    q.push(static_cast<link_type>(x));
    else
    return;
  }
  
  
  
  iterator _insert(base_ptr l, base_ptr r, bool count = true) {
    
#ifdef DEBUG_TREE
    cout<<"inside inserting two for iterators"<<endl;
    cout<<"left pointer "<<l<<endl;
    if (l)
    cout<<"left "<<_value(l)<<endl;
    cout<<"right pointer "<<r<<endl;
    cout<<"right "<<_value(r)<<endl;
#endif
    
    // set parent to header to count nodes
    link_type p = _end();
    r->parent_ = p;
    
    // get left most
    base_ptr lt = r;
    while (lt->left_ != 0)
    lt = lt->left_;
    
#ifdef DEBUG_TREE
    cout<<"left most -> "<<_value(lt)<<endl;
#endif
    
    // count nodes
    if (count)
    for (iterator it(static_cast<link_type>(lt)); it != end(); ++it)
    ++impl_.node_count_;
    
    // get right most
    base_ptr rt = r;
    while (rt->right_ != 0)
    rt = rt->right_;
    
    // if l is not null, create a new node
    link_type lr = static_cast<link_type>(r);
    if (l) {
      
      auto vp = impl_.cost_(_value(l),_value(r));
      lr = _create_node(vp);
      ++impl_.node_count_;
      
      p = _parent(l);
      lr->parent_ = p;
      lr->right_ = r;
      lr->left_ = l;
      
      r->parent_ = lr;
      l->parent_ = lr;
      
    } else
    _leftmost() = lt;
    
    // fix header
    if (lr->parent_ == _end()) {
      
      _root() = lr;
      _rightmost() = rt;
      
    } else {
      // fix parent pointers
      if (p->right_ == l)
      p->right_ = lr;
      else if (p->left_ == l)
      p->left_ = lr;
    }
    
    // recursively fix parents
    auto fix = _value(lr);
    while (lr->parent_ != _end()) {
      p = _parent(lr);
      fix = impl_.cost_(_value(p), fix);
      p->value_ = fix;
      lr = p;
    }
    
    return iterator(static_cast<link_type>(r));
  }
  
  bool _is_leaf(const_base_ptr x) const
  { assert(x != nullptr); return (_right(x) == nullptr) && (_left(x) == nullptr); }
  
  
  
  public:
  
  
  struct _Queue_compare {
    
    // pointer - ancestor expansion - best cost tuples
    typedef std::tuple<link_type, Real, Real> tuple_type;
    
    bool operator()(const tuple_type& t1, const tuple_type& t2) const
    { return std::get<1>(t1) > std::get<1>(t2); }
  };
  
  //! Adapted from Omohundro:89 routine best_sib
  template <typename Arg>
  iterator best_sibling(Arg&& v) {
    
#ifdef DEBUG_TREE
    cout<<"--------------------"<<endl;
    cout<<"inserting "<<v<<endl;
#endif
    
    typedef typename _Queue_compare::tuple_type tuple_type;
    typedef typename std::priority_queue<tuple_type, std::vector<tuple_type>, _Queue_compare> queue_type;
    
    // set result to root
    link_type result = _begin();
    
    // check for empty tree
    if (!result)
    return iterator(nullptr);
    
    // check for leaf node
    else if (_is_leaf(result))
    return iterator(result);
    
    // compute combined cost
    Real bcost = impl_.cost_(impl_.cost_(v, _value(result)));
    
#ifdef DEBUG_TREE
    cout<<"root -> "<<_value(result)<<endl;
    cout<<"bcost -> "<<bcost<<endl;
#endif
    
    // start the queue
    queue_type q;
    decltype(impl_.cost_(v)) e, c;
    
    q.push(std::make_tuple(result, 0., bcost));
    
    // while queue is not empty
    while (!q.empty()) {
      
      // best candidate
      auto t = q.top();
      auto nd = std::get<0>(t);
      auto aexp = std::get<1>(t);
      auto ndvol = std::get<2>(t);
      // remove top most element
      q.pop();
      
#ifdef DEBUG_TREE
      cout<<"  inside while"<<endl;
      cout<<"  top queue item "<<_value(nd)<<endl;
#endif
      
      // no way to get better than bnd
      if (aexp >= bcost) {
#ifdef DEBUG_TREE
        cout<<"*** BREAKING"<<endl;
#endif
        break;
      }
      
      // new ancestor expansion
      e = aexp + ndvol - impl_.cost_(_value(nd));
      
#ifdef DEBUG_TREE
      cout<<"  a exp "<<e<<endl;
      cout<<"  bcost "<<bcost<<endl;
      cout<<"  cost x "<<impl_.cost_(_value(nd))<<endl;
      cout<<"  e -> "<<e<<endl;
#endif
      
      // do left node
      c = impl_.cost_(impl_.cost_(v, _value(nd->left_)));
#ifdef DEBUG_TREE
      cout<<"  c left -> "<<c<<endl;
#endif
      if (c+e <= bcost) {
        bcost = c+e;
        result = _left(nd);
#ifdef DEBUG_TREE
        cout<<"  c+e < bcost!!! result set to left "<<_value(result)<<endl;
#endif
      }
      if(!_is_leaf(nd->left_))
      q.push(std::make_tuple(_left(nd), e, c));
      
      // do right node
      c = impl_.cost_(impl_.cost_(v, _value(nd->right_)));
#ifdef DEBUG_TREE
      cout<<"  c right -> "<<c<<endl;
#endif
      if (c+e <= bcost) {
        bcost = c+e;
        result = _right(nd);
#ifdef DEBUG_TREE
        cout<<"  c+e < bcost!!! result set to right "<<_value(result)<<endl;
#endif
      }
      if(!_is_leaf(nd->right_))
      q.push(std::make_tuple(_right(nd), e, c));
      
    } // while queue is not empty
    
#ifdef DEBUG_TREE
    cout<<"best sibling -> "<<_value(result)<<endl;
#endif
    return iterator(result);
  }
  
  //! Adapted from Omohundro:89 routine cheap_best_sib
  template <typename Arg>
  iterator cheap_best_sibling(Arg&& v) {
    
#ifdef DEBUG_TREE
    cout<<"--------------------"<<endl;
    cout<<"inserting "<<v<<endl;
#endif
    
    link_type x = _begin();
    
    if (!x)
    return iterator(nullptr);
    
#ifdef DEBUG_TREE
    cout<<"root -> "<<_value(x)<<endl;
#endif
    
    link_type result = x;
    
    // compute combined cost
    Real wv = impl_.cost_(v, _value(x));
    // node vol + ancestor expansion
    Real bcost = wv;
    // ancestor expansion starts at zero
    Real ae = 0;
    
#ifdef DEBUG_TREE
    cout<<"bcost -> "<<bcost<<endl;
#endif
    
    // while not a leaf node
    while (x->left_ && x->right_) {
      
#ifdef DEBUG_TREE
      cout<<"  inside while"<<endl;
#endif
      auto nd = impl_.cost_(_value(x));
      
      // correct for both children
      ae += wv - nd;
      
#ifdef DEBUG_TREE
      cout<<"  top queue item "<<_value(x)<<endl;
      cout<<"  a exp "<<ae<<endl;
      cout<<"  bcost "<<bcost<<endl;
      cout<<"  cost x "<<nd<<endl;
      cout<<"  e -> "<<ae<<endl;
#endif
      
      // cannot do any better, insert the volumes in the tree
      if (ae >= bcost)
      break;
      else {
        
        auto lv = impl_.cost_(v, _value(x->left_));
        auto rv = impl_.cost_(v, _value(x->right_));
        
#ifdef DEBUG_TREE
        cout<<"  c left -> "<<lv<<endl;
#endif
        if (ae + lv <= bcost) {
          bcost = ae + lv;
          result = _left(x);
#ifdef DEBUG_TREE
          cout<<"  c+e < bcost!!! result set to left "<<_value(_left(x))<<endl;
#endif
        }
        
#ifdef DEBUG_TREE
        cout<<"  c right -> "<<rv<<endl;
#endif
        
        if (ae + rv <= bcost) {
          bcost = ae + rv;
          result = _right(x);
#ifdef DEBUG_TREE
          cout<<"  c+e < bcost!!! result set to right "<<_value(_right(x))<<endl;
#endif
        }
        
        // if left branch expands less
        if (lv - impl_.cost_(_value(x->left_)) <=
            rv - impl_.cost_(_value(x->right_))) {
          wv = lv;
          x = _left(x);
        } else {
          wv = rv;
          x = _right(x);
        }
      }
    } // while loop
    
#ifdef DEBUG_TREE
    cout<<"best sibling -> "<<_value(result)<<endl;
#endif
    return iterator(result);
  }
  
  private:
  
  template <class functor_type>
  void _execute_at_leaves(iterator x, std::set<base_ptr>& visited, functor_type& fn) {
    
    auto it = visited.find(x.node_);
    if (it != visited.end())
    return;
    
    visited.insert(x.node_);
    
    if (x.is_leaf())
    fn(x);
    else {
      if (x.node_->left_ != 0) _execute_at_leaves(iterator(static_cast<link_type>(x.node_->left_)), visited, fn);
      if (x.node_->right_ != 0) _execute_at_leaves(iterator(static_cast<link_type>(x.node_->right_)), visited, fn);
    }
  }
  
  public:
  
  //! Collect neighbors
  /*! This function is used to collect neighbors from a tree given a functor.
   * The functor should implement a parameterless and a one-parameter operator() overloads.
   * The first load is used as the predicate to finish the search of neighbors, whereas
   * the second one is executed at leaves.
   */
  template <class functor_type>
  void collect_neighbors(iterator it, functor_type& fn) {
    
    // check for null ptr
    if (it.node_ == nullptr)
    return;
    
    std::set<base_ptr> visited;
    
    base_ptr p = it.node_;
    
    // execute while predicate is not true
    while (!fn()) {
      
      p = p->parent_;
      
      // root node
      if (p == _root()->parent_)
      break;
      
      _execute_at_leaves(iterator(static_cast<link_type>(p)), visited, fn);
    }
  }
  
  template <typename Arg>
  std::pair<iterator,bool> cheap_insert(Arg&& v) {
    
    iterator it = cheap_best_sibling(v);
    link_type r = _create_node(v);
    
    return std::pair<iterator, bool>(_insert(it.node_, r), true);
  }
  
  
  template <typename Arg>
  std::pair<iterator,bool> insert(Arg&& v) {
    
    iterator it = best_sibling(v);
    link_type r = _create_node(v);
    
    return std::pair<iterator, bool>(_insert(it.node_, r), true);
  }
  
  
  std::pair<iterator,bool> insert(link_type x)
  { return insert(iterator(x)); }
  
  std::pair<iterator,bool> insert(iterator vit) {
    
    iterator it = best_sibling(_value(vit.node_));
    return std::pair<iterator, bool>(_insert(it.node_, vit.node_), true);
  }
  
  
  iterator sibling(iterator it) {
    
    base_ptr x = it.node_;
    if (!x)
    return iterator(nullptr);
    if (x == _root())
    return end();
    
    base_ptr p = x->parent_;
    assert (p != 0);
    if (p->left_ == x)
    return iterator(_right(p));
    return iterator(_left(p));
  }
  
  public:
  
  iterator relocate_parent(iterator it) {
    
#if DEBUG_TREE
    if (it.node_)
    cout<<"relocating node "<<_value(it.node_)<<endl;
#endif
    
    link_type x = static_cast<link_type>(it.node_);
    
    
    // null pointer
    if (!x)
    return end();
    
    // get parent
    link_type p = _parent(x);
    
    // if parent is header return root
    if (p == _end())
    return iterator(_begin());
    
    // get grand-parent
    link_type gp = _parent(p);
    
    // no need to recompute volume as the parent has
    // the enclosing volume of its children
    
    // sibling
    link_type s = x == p->left_ ? _right(p) : _left(p);
    
#ifdef DEBUG_TREE
    cout<<"relocating x -> "<<_value(x)<<endl;
    cout<<"parent of node -> "<<_value(p)<<endl;
    cout<<"sibling -> "<<_value(s)<<endl;
#endif
    
    if (gp == _end()) {
#ifdef DEBUG_TREE
      cout<<"*** WARNING *** sibling new parent is the header!"<<endl;
#endif
      _root() = s;
    } else if (gp->left_ == p) {
      gp->left_ = s;
#ifdef DEBUG_TREE
      cout<<"setting parent->parent->left to sibling"<<endl;
#endif
    }
    else if (gp->right_ == p) {
      gp->right_ = s;
#ifdef DEBUG_TREE
      cout<<"setting parent->parent->right to sibling"<<endl;
#endif
    } else
    assert(false);
    
    // set parent of sibling
    s->parent_ = p->parent_;
    
    // recursively fix parents
    link_type r = _parent(s);
    while (r != _end()) {
      r->value_ = impl_.cost_(_value(r->right_), _value(r->left_));
      r = _parent(r);
    }
    
    // destroy parent node p and reinsert x into the tree
    x->parent_ = 0;
    p->parent_ = 0;
    p->left_ = 0;
    p->right_ = 0;
    _erase(p);
    
    // get the best sibling
    iterator bit = best_sibling(_value(x));
    
#ifdef DEBUG_TREE
    cout<<"best sibling for relocation -> "<<bit.node_<<endl;
    if (it.node_)
    cout<<" best sibling value -> "<<_value(bit.node_)<<endl;
#endif
    
    iterator result = _insert(bit.node_, x, false);
    
    // set leftmost and rightmost
    base_ptr lt = _root();
    while (lt->left_ != 0)
    lt = lt->left_;
    _leftmost() = lt;
    
    base_ptr rt = _root();
    while (rt->right_ != 0)
    rt = rt->right_;
    _rightmost() = rt;
    
    return result;
  }
  
  
  iterator erase(iterator it) {
    
    link_type x = static_cast<link_type>(it.node_);
    
    // null pointer
    if (!x)
    return end();
    
    // root node
    else if (x == _root()) {
      clear();
      return end();
    }
    
    link_type p = _parent(x);
    
    link_type s; // sibling
    if (x == p->left_) {
      s = _right(p);
      p->right_ = 0; // set to zero to use erase
    } else {
      s = _left(p);
      p->left_ = 0;
    }
    
#ifdef DEBUG_TREE
    cout<<"erasing x -> "<<_value(x)<<endl;
    cout<<"parent of node -> "<<_value(p)<<endl;
    cout<<"sibling -> "<<_value(s)<<endl;
#endif
    
    if (_leftmost() == x) {
      base_ptr l = s;
      while (l->left_ != 0)
      l = l->left_;
      _leftmost() = l;
#ifdef DEBUG_TREE
      cout<<"x is leftmost, so setting leftmost to "<<_value(l)<<endl;
#endif
    } if (_rightmost() == x) {
      base_ptr r = s;
      while (r->right_ != 0)
      r = r->right_;
      _rightmost() = r;
#ifdef DEBUG_TREE
      cout<<"x is rightmost, so setting rightmost to "<<_value(r)<<endl;
#endif
    }
    
    s->parent_ = p->parent_;
    if (p->parent_ == _end()) {
#ifdef DEBUG_TREE
      cout<<"*** WARNING *** sibling new parent is the header!"<<endl;
#endif
      _root() = s;
    } else if (p->parent_->left_ == p) {
      p->parent_->left_ = s;
#ifdef DEBUG_TREE
      cout<<"setting parent->parent->left to sibling"<<endl;
#endif
    }
    else {
      p->parent_->right_ = s;
#ifdef DEBUG_TREE
      cout<<"setting parent->parent->right to sibling"<<endl;
#endif
    }
    
    
#ifdef DEBUG_TREE
    cout<<"p right -> "<<p->right_<<endl;
    cout<<"p left -> "<<p->left_<<endl;
    if (p->right_)
    cout<<"p right value -> "<<_value(p->right_)<<endl;
    if (p->left_)
    cout<<"p left value -> "<<_value(p->left_)<<endl;
#endif
    
    // recursively fix parents
    link_type r = _parent(s);
    while (r != _end()) {
      r->value_ = _value(r->right_) + _value(r->left_);
      r = _parent(r);
    }
    
    // erase subtree rooted at p
    _erase(p);
    
    return iterator(s);
  }
  
  
  
  Real volume() const {
    Real v = 0;
    for (const_iterator it = begin(); it != end(); ++it)
    v += it->measure();
    return v;
  }
  
  
  friend std::ostream& operator<<(std::ostream& os, const Tree& t) {
    
    os<<"Tree: count "<<t.impl_.node_count_<<endl;
    os<<"  header: "<<&t.impl_.header_<<endl;
    if (t.impl_.header_.parent_)
    os<<"  ROOT: "<<t.impl_.header_.parent_<<", value "<<_value(t.impl_.header_.parent_)<<endl;
    os<<"  left: "<<t.impl_.header_.left_<<endl;
    os<<"  right: "<<t.impl_.header_.right_<<endl;
    size_t c = 0;
    Real v = 0;
    os<<endl;
    for (Tree::const_iterator it = t.begin(); it != t.end(); ++it, ++c) {
      
      Real cost = t.impl_.cost_(*it);
      
      if (it.node_ == t.impl_.header_.parent_) {
        os<<"    ROOT "<<it.node_;
        if (it.node_ == t._leftmost())
        os<<" [LEFTMOST]";
        if (it.node_ == t._rightmost())
        os<<" [RIGHTMOST]";
        os<<"\n      value: "<<*it<<endl;
        os<<"      parent: "<<it.node_->parent_<<endl;
      } else {
        os<<"    pointer "<<it.node_;
        if (it.node_ == t._leftmost())
        os<<" [LEFTMOST]";
        if (it.node_ == t._rightmost())
        os<<" [RIGHTMOST]";
        os<<"\n      value: "<<*it<<endl;
        os<<"      parent: "<<it.node_->parent_<<", value "<<_value(it.node_->parent_)<<endl;
      }
      
      if (it.node_->left_)
      os<<"      left: "<<it.node_->left_<<", value "<<_value(it.node_->left_)<<endl;
      if (it.node_->right_)
      os<<"      right: "<<it.node_->right_<<", value "<<_value(it.node_->right_)<<endl;
      if (!it.node_->left_ && !it.node_->right_)
      os<<"      leaf"<<endl;
      os<<"      cost: "<<cost<<endl;
      os<<endl;
      
      v += cost;
    }
    cout<<"total cost -> "<<v<<endl;
    
    if (c != t.impl_.node_count_)
    cout<<"*** WARNING *** Tree node count mismatch "<<t.impl_.node_count_<<" != "<<c<<endl;
    return os;
  }
  
  template <class functor_type>
  functor_type execute_at_level(size_t level, functor_type fn) {
    
    link_type x = _begin();
    if (x == nullptr)
    return fn;
    
    std::queue<link_type> q;
    _add_level(q, x, 0, level);
    
    while (!q.empty()) {
      link_type l = q.front();
      assert(l != 0);
      fn(iterator(l));
      q.pop();
    }
    return fn;
  }
};





template <class object_type>
struct Cost_functor {
  
  auto operator()(const object_type& o) const -> decltype(o.measure())
  { return o.measure(); }
  
  auto operator()(const object_type& o1, const object_type& o2) const -> decltype(o1+o2)
  { return (o1+o2); }
  
};


template <class link_type>
struct Tuple_compare {
  
  typedef std::tuple<link_type, link_type, Real> tuple_type;
  
  bool operator()(const tuple_type& t1, const tuple_type& t2) const
  { return std::get<2>(t1) > std::get<2>(t2); }
};

// bottom up tree construction
template <class tree_type, class volume_container>
tree_type* construct_tree_bottom_up(const volume_container& volumes)
{
  // online tree type definitions
  typedef typename tree_type::value_type volume_type;
  typedef typename tree_type::cost_functor cost_functor;
  typedef typename tree_type::iterator iterator_type;
  
  // queue type definitions
  typedef std::tuple<iterator_type, iterator_type, Real> tuple_type;
  typedef typename std::priority_queue<tuple_type, std::vector<tuple_type>, Tuple_compare<iterator_type> > queue_type;
  
  // already paired container type definitions
  typedef std::set<iterator_type> paired_set;
  typedef typename paired_set::iterator paired_iterator;
  
  size_t numVolumes = volumes.size();
  assert(numVolumes != 0);
  
  tree_type* tp = new tree_type();
  const cost_functor& cost = tp->cost();
  tree_type& t = *tp;
  queue_type q;
  paired_set s;
  
  std::vector<iterator_type> links(numVolumes);
  
  // insert volumes into online tree
  for (size_t i=0; i<numVolumes; ++i) {
    std::pair<iterator_type, bool> p = t.insert(volumes[i]);
    assert(p.second);
    links[i] = p.first;
  }
  
  iterator_type ii,jj;
  for (size_t i=0; i<numVolumes-1; ++i) {
    
    Real merged = std::numeric_limits<Real>::infinity();
    
    for (size_t j=i+1; j<numVolumes; ++j) {
      
      volume_type v = cost(*links[i], *links[j]);
      Real measure = cost(v);
      if (measure < merged) {
        ii = links[i];
        jj = links[j];
        merged = measure;
      }
    }
    q.push(std::make_tuple(ii,jj,merged));
  }
  
  while (!q.empty()) {
    
    // best candidate
    auto tuple = q.top();
    iterator_type node = std::get<0>(tuple);
    iterator_type pair = std::get<1>(tuple);
    
    q.pop();
    
    paired_iterator it = s.find(node);
    if (it != s.end())
    continue;
    
    it = s.find(pair);
    if (it != s.end())
    continue;
    
#if DEBUG_TREE
    cout<<t<<endl;
    cout<<"node -> "<<node.node_<<endl;
    cout<<"node value -> "<<*node<<endl;
    cout<<"pair -> "<<pair.node_<<endl;
    cout<<"pair value -> "<<*pair<<endl;
    cout<<"volume -> "<<std::get<2>(tuple)<<endl;
    
#endif
    
    // obtain best pair from online tree
    iterator_type best = t.sibling(node);
#if DEBUG_TREE
    cout<<"best sibling -> "<<*best<<endl;
#endif
    if (pair == best) {
      
      iterator_type parent = t.relocate_parent(node);
      
      assert(parent != t.end());
      
#if DEBUG_TREE
      cout<<"BEST SIBLING MATCHES! relocating..."<<endl;
      cout<<"  parent -> "<<*parent<<endl;
#endif
      
      // add paired nodes to container
      s.insert(node);
      s.insert(pair);
      
      // compute best pair for the parent
      best = t.sibling(parent);
      
      if (best != t.end()) {
        
        tuple_type tuple = std::make_tuple(parent, best, cost(cost(*parent, *best)));
        q.push(tuple);
        
        iterator_type node = std::get<0>(tuple);
        iterator_type pair = std::get<1>(tuple);
        
#if DEBUG_TREE
        cout<<"  pushing new tuple "<<endl;
        cout<<"  node -> "<<*node<<endl;
        cout<<"  pair -> "<<*pair<<endl;
        cout<<"  volume -> "<<std::get<2>(tuple)<<endl;
#endif
      }
    }
    // else enqueue new pair
    else {
      
      if (best != t.end()) {
        
        tuple_type tuple = std::make_tuple(node,best, cost(cost(*node, *best)));
        q.push(tuple);
        
        iterator_type node = std::get<0>(tuple);
        iterator_type pair = std::get<1>(tuple);
        
#if DEBUG_TREE
        cout<<"BEST SIBLING DOES NOT MATCH!"<<endl;
        cout<<"  pushing new tuple "<<endl;
        cout<<"  node -> "<<*node<<endl;
        cout<<"  pair -> "<<*pair<<endl;
        cout<<"  volume -> "<<std::get<2>(tuple)<<endl;
#endif
      }
      
    }
  }
  
  return tp;
}



template <class volume_type>
class Volume_creator;

template<int d>
struct Volume_creator<Ball<d> > {
  
  typedef typename Ball<d>::point_type point_type;
  
  template <class coord_array>
  static Ball<d> create(const coord_array& coord) {
    
    std::vector<point_type> pts;
    size_t nnodes = coord.size();
    pts.reserve(nnodes);
    for (size_t i=0; i<nnodes; ++i)
    pts.push_back(point_type(coord[i]));
    return bounding_ball<d>(pts);
  }
};

template <int d>
struct Volume_creator<BoundingBox<d> > {
  
  typedef typename BoundingBox<d>::point_type point_type;
  
  template <class coord_array>
  static Sphere create(const coord_array& coord) {
    
    std::vector<point_type> pts;
    size_t nnodes = coord.size();
    pts.reserve(nnodes);
    for (size_t i=0; i<nnodes; ++i)
    pts.push_back(point_type(coord[i]));
    return BoundingBox<d>(pts.begin(), pts.end());
  }
};



// bottom up tree construction
template <class tree_type, class model_type, class element_type>
tree_type* construct_tree_bottom_up(model_type& model) {
  
  // online tree type definitions
  typedef typename tree_type::value_type volume_type;
  typedef typename tree_type::cost_functor cost_functor;
  typedef typename tree_type::iterator iterator_type;
  
  // queue type definitions
  typedef std::tuple<iterator_type, iterator_type, Real> tuple_type;
  typedef typename std::priority_queue<tuple_type, std::vector<tuple_type>, Tuple_compare<iterator_type> > queue_type;
  
  // already paired container type definitions
  typedef std::set<iterator_type> paired_set;
  typedef typename paired_set::iterator paired_iterator;
    
  tree_type* tp = new tree_type();
  tree_type& t = *tp;
  const cost_functor& cost = t.cost();
  std::vector<iterator_type> links;
  
  typedef typename model_type::mesh_type mesh_type;
  
  mesh_type& mesh = model.getMesh();
  int dim = mesh.getSpatialDimension();
  
  // iterate over elements of lower dimension
  typename mesh_type::type_iterator it = mesh.firstType(dim-1);
  typename mesh_type::type_iterator end = mesh.lastType(dim-1);
  
  size_t numVolumes = 0;
  for(; it != end; ++it) {
    
    // add elements to corresponding surface
    for(UInt e = 0; e < mesh.getNbElement(*it); ++e) {
      
      // create solid mechanics element
      element_type el(model, *it, e);
      
      // vector of coordinates for update
      std::vector<const Real*> coord = el.coordinates();
      
      // create leaf volume
      volume_type v = Volume_creator<volume_type>::create(coord);
      
      // do not add zero measure volumes
      if (v.measure() == 0 && dim != 1)
      continue;
      
      std::pair<iterator_type, bool> p = t.insert(v);
      assert(p.second);
      links.push_back(p.first);
      ++numVolumes;
      
      // save leaf data
      if (!t.add_data(p.first,el)) {
        cout<<"*** ERROR *** Could not add data. Aborting..."<<endl;
        exit(1);
      }
      
      
    } // loop over elements
  } // loop over types of element
  
  if (numVolumes == 0) {
    cout<<"*** ERROR *** No volumes were created. Aborting..."<<endl;
    exit(1);
  }
  
  queue_type q;
  paired_set s;
  
  iterator_type ii,jj;
  for (size_t i=0; i<numVolumes-1; ++i) {
    
    Real merged = std::numeric_limits<Real>::infinity();
    
    for (size_t j=i+1; j<numVolumes; ++j) {
      
      volume_type v = cost(*links[i], *links[j]);
      Real measure = v.measure();
      if (measure < merged) {
        ii = links[i];
        jj = links[j];
        merged = measure;
      }
    }
    q.push(std::make_tuple(ii,jj,merged));
  }
  
  while (!q.empty()) {
    
    // best candidate
    auto tuple = q.top();
    iterator_type node = std::get<0>(tuple);
    iterator_type pair = std::get<1>(tuple);
    
    q.pop();
    
    paired_iterator it = s.find(node);
    if (it != s.end())
    continue;
    
    it = s.find(pair);
    if (it != s.end())
    continue;
    
#if DEBUG_TREE
    cout<<t<<endl;
    cout<<"node -> "<<node.node_<<endl;
    cout<<"node value -> "<<*node<<endl;
    cout<<"pair -> "<<pair.node_<<endl;
    cout<<"pair value -> "<<*pair<<endl;
    cout<<"volume -> "<<std::get<2>(tuple)<<endl;
    
#endif
    
    // obtain best pair from online tree
    iterator_type best = t.sibling(node);
#if DEBUG_TREE
    cout<<"best sibling -> "<<*best<<endl;
#endif
    if (pair == best) {
      
      iterator_type parent = t.relocate_parent(node);
      
      assert(parent != t.end());
      
#if DEBUG_TREE
      cout<<"BEST SIBLING MATCHES! relocating..."<<endl;
      cout<<"  parent -> "<<*parent<<endl;
#endif
      
      // add paired nodes to container
      s.insert(node);
      s.insert(pair);
      
      // compute best pair for the parent
      best = t.sibling(parent);
      
      if (best != t.end()) {
        
        tuple_type tuple = std::make_tuple(parent, best, cost(cost(*parent,*best)));
        q.push(tuple);
        
        iterator_type node = std::get<0>(tuple);
        iterator_type pair = std::get<1>(tuple);
        
#if DEBUG_TREE
        cout<<"  pushing new tuple "<<endl;
        cout<<"  node -> "<<*node<<endl;
        cout<<"  pair -> "<<*pair<<endl;
        cout<<"  volume -> "<<std::get<2>(tuple)<<endl;
#endif
      }
    }
    // else enqueue new pair
    else {
      
      if (best != t.end()) {
        
        tuple_type tuple = std::make_tuple(node,best, cost(cost(*node, *best)));
        q.push(tuple);
        
        iterator_type node = std::get<0>(tuple);
        iterator_type pair = std::get<1>(tuple);
        
#if DEBUG_TREE
        cout<<"BEST SIBLING DOES NOT MATCH!"<<endl;
        cout<<"  pushing new tuple "<<endl;
        cout<<"  node -> "<<*node<<endl;
        cout<<"  pair -> "<<*pair<<endl;
        cout<<"  volume -> "<<std::get<2>(tuple)<<endl;
#endif
      }
      
    }
  }
  
  return tp;
}


// bottom up tree construction
template <class tree_type, class element_container>
tree_type* construct_tree_bottom_up(element_container& elems) {
  
  // element type
  typedef typename element_container::value_type element_type;
  
  // online tree type definitions
  typedef typename tree_type::value_type volume_type;
  typedef typename tree_type::cost_functor cost_functor;
  typedef typename tree_type::iterator iterator_type;
  
  // queue type definitions
  typedef std::tuple<iterator_type, iterator_type, Real> tuple_type;
  typedef typename std::priority_queue<tuple_type, std::vector<tuple_type>, Tuple_compare<iterator_type> > queue_type;
  
  // already paired container type definitions
  typedef std::set<iterator_type> paired_set;
  typedef typename paired_set::iterator paired_iterator;
  
  tree_type* tp = new tree_type();
  tree_type& t = *tp;
  const cost_functor& cost = t.cost();
  std::vector<iterator_type> links;
  
  size_t numVolumes = 0;
  
  for (typename element_container::iterator it = elems.begin();
       it != elems.end(); ++it) {
    
    element_type &el = *it;
    
    // vector of coordinates for update
    std::vector<const Real*> coord = el.coordinates();
    
    // create leaf volume
    volume_type v = Volume_creator<volume_type>::create(coord);
    
    // do not add zero measure volumes
    if (v.measure() == 0)
    continue;
    
    std::pair<iterator_type, bool> p = t.insert(v);
    assert(p.second);
    links.push_back(p.first);
    ++numVolumes;
    
    // save leaf data
    if (!t.add_data(p.first,el)) {
      cout<<"*** ERROR *** Could not add data. Aborting..."<<endl;
      exit(1);
    }
    
  } // loop over elements in container
  
  if (numVolumes == 0) {
    cout<<"*** ERROR *** No volumes were created. Aborting..."<<endl;
    exit(1);
  }
  
  queue_type q;
  paired_set s;
  
  iterator_type ii,jj;
  for (size_t i=0; i<numVolumes-1; ++i) {
    
    Real merged = std::numeric_limits<Real>::infinity();
    
    for (size_t j=i+1; j<numVolumes; ++j) {
      
      volume_type v = cost(*links[i], *links[j]);
      Real measure = v.measure();
      if (measure < merged) {
        ii = links[i];
        jj = links[j];
        merged = measure;
      }
    }
    q.push(std::make_tuple(ii,jj,merged));
  }
  
  while (!q.empty()) {
    
    // best candidate
    auto tuple = q.top();
    iterator_type node = std::get<0>(tuple);
    iterator_type pair = std::get<1>(tuple);
    
    q.pop();
    
    paired_iterator it = s.find(node);
    if (it != s.end())
    continue;
    
    it = s.find(pair);
    if (it != s.end())
    continue;
    
#if DEBUG_TREE
    cout<<t<<endl;
    cout<<"node -> "<<node.node_<<endl;
    cout<<"node value -> "<<*node<<endl;
    cout<<"pair -> "<<pair.node_<<endl;
    cout<<"pair value -> "<<*pair<<endl;
    cout<<"volume -> "<<std::get<2>(tuple)<<endl;
    
#endif
    
    // obtain best pair from online tree
    iterator_type best = t.sibling(node);
#if DEBUG_TREE
    cout<<"best sibling -> "<<*best<<endl;
#endif
    if (pair == best) {
      
      iterator_type parent = t.relocate_parent(node);
      
      assert(parent != t.end());
      
#if DEBUG_TREE
      cout<<"BEST SIBLING MATCHES! relocating..."<<endl;
      cout<<"  parent -> "<<*parent<<endl;
#endif
      
      // add paired nodes to container
      s.insert(node);
      s.insert(pair);
      
      // compute best pair for the parent
      best = t.sibling(parent);
      
      if (best != t.end()) {
        
        tuple_type tuple = std::make_tuple(parent, best, cost(cost(*parent,*best)));
        q.push(tuple);
        
        iterator_type node = std::get<0>(tuple);
        iterator_type pair = std::get<1>(tuple);
        
#if DEBUG_TREE
        cout<<"  pushing new tuple "<<endl;
        cout<<"  node -> "<<*node<<endl;
        cout<<"  pair -> "<<*pair<<endl;
        cout<<"  volume -> "<<std::get<2>(tuple)<<endl;
#endif
      }
    }
    // else enqueue new pair
    else {
      
      if (best != t.end()) {
        
        tuple_type tuple = std::make_tuple(node,best, cost(cost(*node, *best)));
        q.push(tuple);
        
        iterator_type node = std::get<0>(tuple);
        iterator_type pair = std::get<1>(tuple);
        
#if DEBUG_TREE
        cout<<"BEST SIBLING DOES NOT MATCH!"<<endl;
        cout<<"  pushing new tuple "<<endl;
        cout<<"  node -> "<<*node<<endl;
        cout<<"  pair -> "<<*pair<<endl;
        cout<<"  volume -> "<<std::get<2>(tuple)<<endl;
#endif
      }
      
    }
  }
  
  return tp;
}




template <class tree_type>
void print_mathematica(tree_type& t) {
  
  typedef typename tree_type::value_type volume_type;
  
  int d = volume_type::dim();
  
  cout<<std::fixed;
  cout.precision(2);
  
  if (d == 2)
  cout<<"Graphics[{";
  else
  cout<<"Graphics3D[{";
  
  for (typename tree_type::const_leaf_iterator it = t.leaves_begin(); it != t.leaves_end(); ++it)
  cout<<*it<<", ";
  
  cout<<"}];"<<endl;
}


template <typename T>
std::string stringify(const T& x) {
  std::stringstream o;
  if (!(o << x)) {
    cout<<"*** ERROR *** Bad argument to stringity function"<<endl;
    exit(1);
  }
  return o.str();
}

template <class tree>
void export_video(tree& t, Real duration, Real min_opacity, std::string pre = "", std::string post = "") {
  
  typedef typename tree::value_type volume_type;
  typedef typename volume_type::aabb_type aabb_type;
  typedef typename aabb_type::point_type point_type;
  typedef typename point_type::value_type value_type;
  typedef typename tree::iterator iterator;
  
  int d = volume_type::dim();
  cout<<std::fixed;
  cout.precision(2);
  
  std::tuple<size_t, size_t, size_t> h = t.height();
  
  typename tree::iterator it = t.begin();
  
  aabb_type bb = it->bounding_box();
  
  // get bounding box
  for (; it != t.end(); ++it)
  bb += it->bounding_box();
  
  std::string figs;
  int H = std::get<2>(h);
  int k = 0;
  for (int l = H; l>=0; --l) {
    
    std::string g = "g";
    g += stringify(k);
    ++k;
    
    // compute opacity for figure
    Real opacity = (1 - min_opacity)*static_cast<double>(l)/H + min_opacity;
    
    if (d == 2)
    cout<<g<<" = Graphics[{Opacity["<<stringify(opacity)<<"],"<<pre;
    else
    cout<<g<<" = Graphics3D[{Opacity["<<stringify(opacity)<<"],"<<pre;
    
    t.execute_at_level(l, [](iterator it){ cout<<*it<<", ";});
    
    figs += figs.empty() ? g : (','+g);
    
    cout<<"}];\nShow["<<figs<<", PlotRange-> {{";
    const point_type &m = bb.min();
    const point_type &M = bb.max();
    for (int i=0; i<d; ++i) {
      value_type l = 0.05*(M[i] - m[i]);
      cout<<(m[i]-l)<<','<<(M[i]+l)<<(i < d-1 ? "},{" : "}}");
    }
    cout<<post<<"];\n";
    cout<<"Export[\""<<g<<".png\", %, ImageResolution -> 300];"<<endl;
  }
  
  cout<<"files =  Last /@ Sort[{Characters@#, #} & /@ FileNames[\"*.png\"]]"<<endl;
  cout<<"images = Map[Import, files];"<<endl;
  cout<<"Export[\"video.mov\", images, \"FrameRate\" -> "<<stringify(static_cast<double>(k)/duration)<<"];"<<endl;
}


//// bottom up tree construction
//template <class volume_container>
//Tree<typename volume_container::value_type, Cost_functor<typename volume_container::value_type> >
//construct_tree_bottom_up2(const volume_container& volumes)
//{
//
//  typedef typename volume_container::value_type volume_type;
//
//  // online tree type definitions
//  typedef Tree<volume_type, Cost_functor<volume_type> > tree_type;
//
//  typedef Tree_node<volume_type> node_type;
//  typedef typename tree_type::link_type link_type;
//
//  size_t numVolumes = volumes.size();
//  assert(numVolumes != 0);
//
//  // allocate leafs of the tree
//  link_type* leafs = new link_type[numVolumes];
//
//  for (size_t i=0; i<numVolumes; ++i)
//    leafs[i] = new node_type(volumes[i]);
//
//  while (numVolumes > 1) {
//
//    // find indices of the two "nearest" nodes, based on some criterion
//    size_t ii,jj;
//    Real merged = std::numeric_limits<Real>::infinity();
//    for (size_t i=0; i<numVolumes; ++i)
//      for (size_t j=i+1; j<numVolumes; ++j) {
//
//        volume_type v = leafs[i]->value_ + leafs[j]->value_;
//        if (v.measure() < merged) {
//          ii = i;
//          jj = j;
//          merged = v.measure();
//        }
//      }
//
//    volume_type v = leafs[ii]->value_ + leafs[jj]->value_;
//    link_type parent = new node_type(v);
//
//#if DEBUG_TREE
//    cout<<"best volume found -> "<<v<<endl;
//    cout<<"measure -> "<<v.measure()<<endl;
//    cout<<"parent -> "<<parent->value_<<endl;
//    cout<<"left -> "<<leafs[ii]->value_<<endl;
//    cout<<"right -> "<<leafs[jj]->value_<<endl;
//#endif
//
//    parent->left_ = leafs[ii];
//    parent->right_ = leafs[jj];
//    leafs[ii]->parent_ = parent;
//    leafs[jj]->parent_ = parent;
//
//    // remove the two nodes from the active set and add in the new node.
//    leafs[ii] = parent;
//    leafs[jj] = leafs[numVolumes - 1];
//    --numVolumes;
//  }
//
//  link_type root = leafs[0];
//
//  tree_type t(root);
//
//  delete leafs;
//
//  return t;
//}

__END_AKANTU__

#endif /* __AKANTU_TREE_HH__ */
