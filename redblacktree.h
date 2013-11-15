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

#ifndef INTERVALTREE_H
#define INTERVALTREE_H

#include <iostream>
#include <cstdlib>

namespace torali
{


  class IntervalTreeNode {
  public:
    int key;
    bool color;
    IntervalTreeNode* left;
    IntervalTreeNode* right;
    IntervalTreeNode* p;

    IntervalTreeNode() : key(0), color(false), left(NULL), right(NULL), p(NULL) {}

    IntervalTreeNode(int k) : key(k), color(false), left(NULL), right(NULL), p(NULL) {}
  };

  class IntervalTree {
  public:
    IntervalTreeNode* nil;
    IntervalTreeNode* root;
    
    IntervalTree() {
      nil = new IntervalTreeNode();
      root = nil;
    }

    ~IntervalTree() {
      delete nil;
    }

    void insertFixup(IntervalTreeNode* z) {
      IntervalTreeNode* y;
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
	

    void insertNode(IntervalTreeNode* z) {
      IntervalTreeNode* y = nil;
      IntervalTreeNode* x = root;
      while (x != nil) {
	y = x;
	if (z->key < x->key) x = x->left;
	else x = x->right;
      }
      z->p = y;
      if (y == nil) root = z;
      else {
	if (z->key < y->key) y->left = z;
	else y->right = z;
      }
      z->left = nil;
      z->right = nil;
      z->color = true;
      insertFixup(z);
    }

    void insertKey(int k) {
      insertNode(new IntervalTreeNode(k));
    }

    IntervalTreeNode* searchNode(int k) {
      IntervalTreeNode* x = root;
      while ((x != nil) && (k != x->key)) {
	if (k < x->key) x = x->left;
	else x = x->right;
      }
      return x;
    }

    IntervalTreeNode* minimumNode(IntervalTreeNode* x) {
      while(x->left != nil) x = x->left;
      return x;
    }

    IntervalTreeNode* maximumNode(IntervalTreeNode* x) {
      while(x->right != nil) x = x->right;
      return x;
    }

    IntervalTreeNode* successorNode(IntervalTreeNode* x) {
      if (x->right != nil) return minimumNode(x->right);
      IntervalTreeNode* y =  x->p;
      while ((y != nil) && (x == y->right)) {
	x = y;
	y = y->p;
      }
      return y;
    }

    void leftRotate(IntervalTreeNode* x) {
      IntervalTreeNode* y = x->right;
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
    }

    void rightRotate(IntervalTreeNode* y) {
      IntervalTreeNode* x = y->left;
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
    }
    

    void deleteFixup(IntervalTreeNode* x) {
      IntervalTreeNode* w;
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

    void deleteNode(IntervalTreeNode* z) {
      IntervalTreeNode* y;
      if ((z->left == nil) || (z->right == nil)) y = z;
      else y = successorNode(z);
      IntervalTreeNode* x;
      if (y->left != nil) x = y->left;
      else x = y->right;
      x->p = y->p;
      if (y->p == nil) root = x;
      else {
	if (y == y->p->left) y->p->left = x;
	else y->p->right = x;
      }
      if (y != z) z->key = y->key;
      if (!y->color) deleteFixup(x);
      delete y;
    }

    void deleteKey(int k) {
      deleteNode(searchNode(k));
    }


    void inorder(IntervalTreeNode* x) {
      if (x != nil) {
	inorder(x->left);
	std::cout << x->key << " (P: " << x->p->key << ", L: " << x->left->key << ", R: " << x->right->key << ", C: " << x->color << "); ";
	inorder(x->right);
      }
    }

  };

}

#endif
