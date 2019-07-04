/*
Author: Agnes Valenti
This file contains two functions to obtain adjacent vertices and plaquettes, respectevily to a given spin
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <random>
#include <iomanip>
#include <limits>
#include <ctime>
using namespace std;

void get_2vertices(int spin, int & s1, int & s2, int L){
   //returns two vertices s1 and s2 adjacent to spin
   int spin_x,spin_y,spin_AB;
   spin_x=spin/(2*L);
   spin_y=(spin-spin_x*2*L)/2;
   spin_AB=(spin-spin_x*2*L-spin_y*2);
   int xs1,xs2,ys1,ys2;
   if (spin_AB==0){
      ys1=spin_y;
      ys2=spin_y;
      if (spin_x==0)
        xs1=L-1;
      else
        xs1=(spin_x-1);
      xs2=spin_x;
      }
   if (spin_AB==1){
      xs1=spin_x;
      xs2=spin_x;
      ys1=spin_y;
      ys2=(spin_y+1)%L;
      }
   s1=L*xs1+ys1;
   s2=L*xs2+ys2;
   }

inline void get_2plaquettes(int spin, int & p1, int & p2, int L){
   //returns two plaquettes p1 and p2 adjacent to spin
   int spin_x,spin_y,spin_AB;
   spin_x=spin/(2*L);
   spin_y=(spin-spin_x*2*L)/2;
   spin_AB=(spin-spin_x*2*L-spin_y*2);
   int xp1,xp2,yp1,yp2;
   if (spin_AB==1){
      yp1=spin_y;
      yp2=spin_y;
      xp1=spin_x;
      xp2=(spin_x+1)%L;
      }
   if (spin_AB==0){
      xp1=spin_x;
      xp2=spin_x;
      yp1=spin_y;
      if (spin_y==0)
        yp1=L-1;
      else
        yp1=(spin_y-1);
      yp2=spin_y%L;
      }
   p1=L*xp1+yp1;
   p2=L*xp2+yp2;
   }
