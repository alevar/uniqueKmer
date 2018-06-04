#include <cmath>
#include <string>
#include <iostream>

#include "CM.h"

CM::CM(){

}

CM::CM(float epsilon,float delta){
	this->width=ceil(exp(1)/epsilon);
	this->depth=ceil(log(1/delta));

	sketch = new uint8_t * [depth];
    for (size_t i = 0; i < this->depth; i++) {
        sketch[i] = new uint8_t [width];
        for (size_t j=0;j<this->width;j++){
        	sketch[i][j]=0;
        }
    }

    // std::cout<<"width: "<<width<<" depth: "<<depth<<std::endl;
}

CM::~CM(){
	for (size_t i=0;i<this->depth;i++) {
        delete[] this->sketch[i];
    }
    delete[] this->sketch;
}

void CM::set_params(float epsilon,float delta){
	this->width=ceil(exp(1)/epsilon);
	this->depth=ceil(log(1/delta));

	sketch = new uint8_t * [depth];
    for (size_t i = 0; i < this->depth; i++) {
        sketch[i] = new uint8_t [width];
        for (size_t j=0;j<this->width;j++){
        	sketch[i][j]=0;
        }
    }

    // std::cout<<"width: "<<width<<" depth: "<<depth<<std::endl;
}

void CM::update(std::string& item){
	for (size_t i=0;i<this->depth;i++){
		this->hash=this->hasher(item+std::to_string(i))%this->width;
		this->curVal=sketch[i][this->hash];
		if (curVal<MAX_PERMITTED_COUNT){
			sketch[i][this->hash]=curVal+1;
		}
	}
}

uint8_t CM::estimate(std::string& item){
	uint8_t minCount=MAX_PERMITTED_COUNT;

	for (size_t i=0;i<this->depth;i++){
		this->hash=this->hasher(item+std::to_string(i))%this->width;
		this->curVal=sketch[i][this->hash];
		if (this->curVal<minCount){
			minCount=this->curVal;
		}
	}
	return minCount;
}

// merge another sketch into current
void CM::merge(CM* cm){
	for (size_t i=0;i<this->depth;i++){
		for (size_t j=0;j<this->width;j++){
			// std::cout<<"val at: ["<<j<<","<<i<<"]"<<std::endl;
			this->curVal=this->sketch[i][j]+cm->sketch[i][j];
			if (this->sketch[i][j]<MAX_PERMITTED_COUNT-cm->sketch[i][j]){
				this->sketch[i][j]+=((cm->sketch)[i][j]);
			}	
		}
	}
}

void CM::save(FILE * fp){
	std::cout<<std::endl;
}

void CM::load(FILE * fp){
	std::cout<<std::endl;
}
