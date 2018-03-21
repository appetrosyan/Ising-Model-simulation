//
//  main.cpp
//  Ising Model
//
//  Created by Aleksandr Petrosyan on 20/03/2018.
//  Copyright © 2018 Aleksandr Petrosyan. All rights reserved.
//

#include <iostream>

using namespace std;

class lattice{
    int size;
    short* spin;
    double J;
    double H;
public:
    lattice(int _size, double _J, double _H){
        this->size = _size;
        this->spin = new short(size*size);
        for(int i=0; i<size*size; i++){
            spin[i] = 1;
        }
        this->J = _J;
        this->H = _H;
    }
    
    void print(){
        cout<<size<<endl<<endl;
        int area = size*size;
        for(int i=0; i< area; i++){
            cout<<spin[i]<<' ';
            if (i%size ==size-1 ) cout<<endl;
        }
        cout<<endl;
    }
    
    short get(int row, int col){
        if(row<0 || row >=size){
            cerr<<"Out of bounds"<<row;
            return 0;
        }
        if(col<0 || col>=size){
            cerr<<"Out of bounds"<<col;
            return 0;
        }
        return spin[row*size+col];
    }
    
    void flip(int row, int col){
        if(row<0 || row >=size){
            cerr<<"Out of bounds"<<row;
            return;
        }
        if(col<0 || col>=size){
            cerr<<"Out of bounds"<<col;
            return;
        }
        
    }
    
    int get_size(){
        return size;
    }
    
    float compute_point_energy(int row, int col){
        int accumulator=0;
        if (row >= size || col >= size){
            fprintf(stderr, "Out of bounds");
            return 0.0;
        }
        
        if(row > 0){
            accumulator+=get(row-1, col);
            cerr<<"check"<<accumulator<<endl;
        } else{
            accumulator+=get(size-1,col);
        }
        if (row < size - 1){
            accumulator+=get(row+1, col);
            cerr<<"check"<<accumulator<<endl;
        } else {
            accumulator+=get(0, col);
        }
        
        if (col<0){
            accumulator+=get(row, col-1);
            cerr<<"check"<<accumulator<<endl;
        } else {
            accumulator+=get(row, size-1);
        }
        if (col< size-1) {
            accumulator+=get(row, col+1);
            cerr<<"check"<<accumulator<<endl;
        } else {
            accumulator+=get(row, 0);
        }
        cerr<<accumulator<<" "<<J<<endl;
        return accumulator*J;
    }
    
};

int main(int argc, const char * argv[]) {
    lattice l = lattice (4, 0.5, 0);
    l.print();
    cout<<endl;
    cout<<l.compute_point_energy(0, 0)<<endl;
    return 0;
}

