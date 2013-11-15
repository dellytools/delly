/*
============================================================================
DELLY: Structural variant discovery by integrated PE mapping and SR analysis
============================================================================
Copyright (C) 2012 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/


// Intervals need to be unique!!!

#ifndef INTERVALTREE_H
#define INTERVALTREE_H

#include <iostream>
#include <cstdlib>
#include <vector>

namespace torali
{

  template<typename TCargo>
  class Interval {
  public:
    int low;
    int high;
    TCargo cargo;

    Interval() : low(0), high(0) {}
    Interval(int a, int b) : low(a), high(b) {}
    Interval(int a, int b, TCargo& c) : low(a), high(b), cargo(c) {}
    Interval(int a, int b, TCargo const& c) : low(a), high(b), cargo(c) {}

    Interval(Interval& i) {
      low = i.low;
      high = i.high;
      cargo = i.cargo;
    }
    Interval(Interval const& i) {
      low = i.low;
      high = i.high;
      cargo = i.cargo;
    }
	
  };

  template<typename TInterval>
  class IntervalTreeNode {
  public:
    int maximum;
    bool color;
    IntervalTreeNode* left;
    IntervalTreeNode* right;
    IntervalTreeNode* p;
    TInterval interv;

    IntervalTreeNode() : maximum(0), color(false), left(NULL), right(NULL), p(NULL) {}
    IntervalTreeNode(TInterval& i) : maximum(0), color(false), left(NULL), right(NULL), p(NULL), interv(i) {}
    IntervalTreeNode(TInterval const& i) : maximum(0), color(false), left(NULL), right(NULL), p(NULL), interv(i) {}
  };

  template<typename TInterval>
  class IntervalTree {
  public:
    IntervalTreeNode<TInterval>* nil;
    IntervalTreeNode<TInterval>* root;
    
    IntervalTree() {
      nil = new IntervalTreeNode<TInterval>();
      root = nil;
    }

    ~IntervalTree() {
      while (root != nil) deleteNode(root);
      delete nil;
    }

    void insertFixup(IntervalTreeNode<TInterval>* z) {
      IntervalTreeNode<TInterval>* y;
      while (z->p->color) {
	if (z->p == z->p->p->left) {
	  y = z->p->p->right;
	  if (y->color) {
	    z->p->color = false;
	    y->color = false;
	    z->p->p->color = true;
	    z = z->p->p;
	  } else {
	    if (z == z->p->right) {
	      z = z->p;
	      leftRotate(z);
	    }
	    z->p->color = false;
	    z->p->p->color = true;
	    rightRotate(z->p->p);
	  }
	} else {
	  y = z->p->p->left;
	  if (y->color) {
	    z->p->color = false;
	    y->color = false;
	    z->p->p->color = true;
	    z = z->p->p;
	  } else {
	    if (z == z->p->left) {
	      z = z->p;
	      rightRotate(z);
	    }
	    z->p->color = false;
	    z->p->p->color = true;
	    leftRotate(z->p->p);
	  }
	}
      }
      root->color = false;
    }    
	

    void insertNode(IntervalTreeNode<TInterval>* z) {
      IntervalTreeNode<TInterval>* y = nil;
      IntervalTreeNode<TInterval>* x = root;
      while (x != nil) {
	y = x;
	if ( ((z->interv.low != x->interv.low) ? (z->interv.low < x->interv.low) : (z->interv.high < x->interv.high) ) ) x = x->left;
	else x = x->right;
      }
      z->p = y;
      if (y == nil) root = z;
      else {
	if ( ((z->interv.low != y->interv.low) ? (z->interv.low < y->interv.low) : (z->interv.high < y->interv.high) ) ) y->left = z;
	else y->right = z;
      }
      z->left = nil;
      z->right = nil;
      z->color = true;
      z->maximum = z->interv.high;
      x = z->p;
      while (x != nil) {
	x->maximum = std::max(x->interv.high , std::max(x->left->maximum, x->right->maximum));
	x = x->p;
      }
      insertFixup(z);
    }

    void insertInterval(TInterval& i) {
      insertNode(new IntervalTreeNode<TInterval>(i));
    }

    IntervalTreeNode<TInterval>* searchIntervalExact(TInterval& i) {
      IntervalTreeNode<TInterval>* x = root;
      while ((x != nil) && ((i.low != x->interv.low) || (i.high != x->interv.high))) {
	if ( ((i.low != x->interv.low) ? (i.low < x->interv.low) : (i.high < x->interv.high) ) ) x = x->left;
	else x = x->right;
      }
      return x;
    }

    void enumOverlapInterval(IntervalTreeNode<TInterval>* x, TInterval& i, std::vector<TInterval>& result) {
      if ((x != nil) && (i.low <= x->maximum)) {
	enumOverlapInterval(x->left, i, result);
	if ((i.low <= x->interv.high) && (i.high >= x->interv.low)) result.push_back(x->interv);
	if (i.high >= x->interv.low) enumOverlapInterval(x->right, i, result);
      }
    }	  

    void enumOverlapInterval(TInterval& i, std::vector<TInterval>& result) {
      enumOverlapInterval(root, i, result);
    }

    IntervalTreeNode<TInterval>* minimumNode(IntervalTreeNode<TInterval>* x) {
      while(x->left != nil) x = x->left;
      return x;
    }

    IntervalTreeNode<TInterval>* maximumNode(IntervalTreeNode<TInterval>* x) {
      while(x->right != nil) x = x->right;
      return x;
    }

    IntervalTreeNode<TInterval>* successorNode(IntervalTreeNode<TInterval>* x) {
      if (x->right != nil) return minimumNode(x->right);
      IntervalTreeNode<TInterval>* y =  x->p;
      while ((y != nil) && (x == y->right)) {
	x = y;
	y = y->p;
      }
      return y;
    }

    void leftRotate(IntervalTreeNode<TInterval>* x) {
      IntervalTreeNode<TInterval>* y = x->right;
      x->right = y->left;
      if (y->left != nil) y->left->p = x;
      y->p = x->p;
      if (x->p == nil) root = y;
      else {
	if (x == x->p->left) x->p->left = y;
	else x->p->right = y;
      }
      y->left = x;
      x->p = y;
      x->maximum = std::max(x->interv.high , std::max(x->left->maximum, x->right->maximum));
      y->maximum = std::max(y->interv.high , std::max(y->left->maximum, y->right->maximum));
    }

    void rightRotate(IntervalTreeNode<TInterval>* y) {
      IntervalTreeNode<TInterval>* x = y->left;
      y->left = x->right;
      if (x->right != nil) x->right->p = y;
      x->p = y->p;
      if (y->p == nil) root = x;
      else {
	if (y == y->p->left) y->p->left = x;
	else y->p->right = x;
      }
      x->right = y;
      y->p = x;
      y->maximum = std::max(y->interv.high , std::max(y->left->maximum, y->right->maximum));
      x->maximum = std::max(x->interv.high , std::max(x->left->maximum, x->right->maximum));
    }


    void deleteFixup(IntervalTreeNode<TInterval>* x) {
      IntervalTreeNode<TInterval>* w;
      while ((x != root) && (!x->color)) {
	if (x == x->p->left) {
	  w = x->p->right;
	  if (w->color) {
	    w->color = false;
	    x->p->color = true;
	    leftRotate(x->p);
	    w = x->p->right;
	  }
	  if ((!w->left->color) && (!w->right->color)) {
	    w->color = true;
	    x = x->p;
	  } else {
	    if (!w->right->color) {
	      w->left->color = false;
	      w->color = true;
	      rightRotate(w);
	      w = x->p->right;
	    }
	    w->color = x->p->color;
	    x->p->color = false;
	    w->right->color = false;
	    leftRotate(x->p);
	    x = root;
	  }
	} else {
	  w = x->p->left;
	  if (w->color) {
	    w->color = false;
	    x->p->color = true;
	    rightRotate(x->p);
	    w = x->p->left;
	  }
	  if ((!w->right->color) && (!w->left->color)) {
	    w->color = true;
	    x = x->p;
	  } else {
	    if (!w->left->color) {
	      w->right->color = false;
	      w->color = true;
	      leftRotate(w);
	      w = x->p->left;
	    }
	    w->color = x->p->color;
	    x->p->color = false;
	    w->left->color = false;
	    rightRotate(x->p);
	    x = root;
	  }
	}
      }
      x->color = false;
    }

    void deleteNode(IntervalTreeNode<TInterval>* z) {
      IntervalTreeNode<TInterval>* y;
      if ((z->left == nil) || (z->right == nil)) y = z;
      else y = successorNode(z);
      IntervalTreeNode<TInterval>* x;
      if (y->left != nil) x = y->left;
      else x = y->right;
      x->p = y->p;
      if (y->p == nil) root = x;
      else {
	if (y == y->p->left) y->p->left = x;
	else y->p->right = x;
      }
      if (y != z) z->interv = y->interv;
      if ( x != nil) x->maximum = std::max(x->interv.high , std::max(x->left->maximum, x->right->maximum)); 
      IntervalTreeNode<TInterval>* walker = x->p;
      while (walker != nil) {
	walker->maximum = std::max(walker->interv.high , std::max(walker->left->maximum, walker->right->maximum));
	walker = walker->p;
      }
      if (!y->color) deleteFixup(x);
      delete y;
    }

    void deleteInterval(TInterval& i) {
      deleteNode(searchIntervalExact(i));
    }

    void inorder(IntervalTreeNode<TInterval>* x) {
      if (x != nil) {
	inorder(x->left);
	std::cout << "[" << x->interv.low << "," << x->interv.high << "] (P: " << x->p->interv.low << "," << x->p->interv.high << "; L: " << x->left->interv.low << "," << x->left->interv.high << "; R: " << x->right->interv.low << "," << x->right->interv.high << "; C: " << x->color << "; M: " << x->maximum << "); ";
	//if (x->maximum != std::max(x->interv.high , std::max(x->left->maximum, x->right->maximum))) {
	//std::cout << "Maximum wrong!" << std::endl;
	//exit(-1);
	//}
	inorder(x->right);
      }
    }

    void inorder() {
      inorder(root);
    }

  };

}

#endif
