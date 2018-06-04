#ifndef INCL_CM
#define INCL_CM

#include <string>

#define MAX_PERMITTED_COUNT 65535

class CM {

	private:
		unsigned short int** sketch,curVal;
		unsigned int width,depth;
		std::hash<std::string> hasher;
		size_t hash;

	public:
      CM();
  		CM(float epsilon,float delta);
  		~CM();
      
      void set_params(float epsilon,float delta);
  		void update(std::string& item);
  		void merge(CM * cm);
  		void save(FILE * fp);
  		void load(FILE * fp);
  		unsigned short int estimate(std::string& item);
};

#endif