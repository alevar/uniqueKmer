#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <algorithm>
#include <dirent.h>
#include <stdlib.h>
#include <vector>
#include <omp.h>

#include "CM.h"
#include "arg_parse.h"

// ./uniqueLoci_build -i ./testData -o ./log31.csv -k 31 -e 0.00000001 -d 0.005 -t 2

typedef std::vector< std::pair<long,std::string> > fsMap; //type of object that holds file paths and sizes
typedef std::vector< fsMap > thread_fsMap; //type of object that holds file paths and sizes for each thread

std::string reverse(std::string seq){ //https://stackoverflow.com/questions/33074574/creating-complement-of-dna-sequence-and-reversing-it-c
    bool ambiguous=false;
    auto lambda = [&ambiguous](const char c) {
        switch (c) {
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        case 'N':
            return 'N';
        default:
            // std::cout<<"ERROR: "<<c<<std::endl;
            // throw std::domain_error("Invalid nucleotide.");
            ambiguous=true;
            return 'N';
        }
    };

    std::transform(seq.cbegin(), seq.cend(), seq.begin(), lambda);
    if (ambiguous==false){
        return seq;
    }
    return {};
}

// UNALIGN SAM'S code
enum Opt {INPUT    = 'i',
          OUTPUT       = 'o',
          KMER_LEN     = 'k',
          EPSILON      = 'e',
          DELTA        = 'd',
          THREADS      = 't',
          STEP         = 's',
          MAX_COUNT    = 'c'};

CM combineCM(CM r,CM n) {
    // r is the already reduced CM
    // n is the new CM
    r.merge(&n);
    return r;
}

int main(int argc, char** argv) {

    ArgParse args("Find unique loci in a genome");

    args.add_string(Opt::INPUT, "ref_fasta", "", "directory where reference sequence directories are stored.");
    args.add_string(Opt::OUTPUT, "output", "", "directory where the database will be stored");
    args.add_int(Opt::KMER_LEN, "kmer_length", 31, "kmer length to use");
    args.add_double(Opt::EPSILON, "epsilon", 0.00000001, "epsilon parameter of the CM sketch");
    args.add_double(Opt::DELTA, "delta", 0.001, "delta parameter of CM sketch");
    args.add_int(Opt::THREADS, "threads",1,"number of threads to run concurrently");
    args.add_int(Opt::STEP, "step",1,"how many positions to skip");
    args.add_int(Opt::MAX_COUNT, "max_count",1,"Only kmers with a count <= to max_ount will be reported in the log");
    args.parse_args(argc, argv);

    std::string inputDir=args.get_string(Opt::INPUT);

    int kmerLen=args.get_int(Opt::KMER_LEN);
    double epsilon=args.get_double(Opt::EPSILON);
    double delta=args.get_double(Opt::DELTA);
    unsigned int numThreads=args.get_int(Opt::THREADS);
    int max_count=args.get_int(Opt::MAX_COUNT);

    /*================================================================================
    =========================FIND FILES AND ASSIGN TO THREADS=========================
    ================================================================================*/
    DIR *parentDir=opendir(inputDir.c_str());
    struct dirent *parentEnt;
    std::string parentFP=inputDir;
    std::string subFP=inputDir;
    if (strcmp(parentFP.substr(parentFP.length()-1).c_str(),"/")==0){
        parentFP=parentFP.substr(0,parentFP.length()-1);
    }

    std::string curFP;
    fsMap fsm;

    if (parentDir!=NULL){
        while ((parentEnt = readdir (parentDir)) != NULL) {
            if (strcmp(parentEnt->d_name,".")!=0 && strcmp(parentEnt->d_name,"..")!=0){
                subFP=parentFP;
                subFP+="/";
                subFP+=parentEnt->d_name;
                DIR *subDir=opendir(subFP.c_str());
                struct dirent *subEnt;
                if (subDir!=NULL){
                    while ((subEnt = readdir (subDir)) != NULL) {
                        if (strcmp(subEnt->d_name,".")!=0 && strcmp(subEnt->d_name,"..")!=0 && std::string(subEnt->d_name).length()>4){
                            curFP=std::string(subEnt->d_name);
                            if (strcmp(curFP.substr(curFP.length()-4).c_str(),".fna")==0 || strcmp(curFP.substr(curFP.length()-3).c_str(),".fa")==0){
                                std::ifstream* in=new std::ifstream((subFP+"/"+curFP).c_str(), std::ifstream::ate | std::ifstream::binary);
                                fsm.push_back(std::pair<long,std::string>(in->tellg(),subFP+"/"+curFP));
                                delete in;
                            }
                        }
                    }
                }
                closedir (subDir);
            }
        }
    }
    else{
        std::cout<<"Error opening parent directory"<<std::endl;
        exit (EXIT_FAILURE);
    }

    // if threads requested is greater than the number of files - reduce number of threads to the number of files
    if (numThreads>fsm.size()){
        numThreads=fsm.size();
    }

    // initialize thread assignments
    thread_fsMap threadFPs;
    std::sort(fsm.begin(), fsm.end());

    // first distribute n largest objects into threads
    for (unsigned int i=0;i<numThreads;i++){
        threadFPs.push_back(fsMap());
        threadFPs[i].push_back(fsm.back());
        fsm.pop_back();
    }

    // next begin distributing smallest objects to the smallest bins
    // starting at the smallest bin
    if (fsm.size()>0){
        while(true){
            std::sort(threadFPs.begin(),threadFPs.end(),[](fsMap a, fsMap b){
                long s1=0;
                long s2=0;
                for (std::pair<long,std::string> i : a)
                    s1+=i.first;
                for (std::pair<long,std::string> i : b)
                    s2+=i.first;
                return s1 < s2;
            });
            threadFPs[0].push_back(fsm.front());
            fsm.erase(fsm.begin());
            if (fsm.size()==0){
                break;
            }
        }
    }
    closedir (parentDir);

    // for (int i=0;i<threadFPs.size();i++){
    //     std::cout<<"\nTHREAD #"<<i<<std::endl;
    //     for (auto j : threadFPs[i])
    //         std::cout<<"File: "<<j.second<<" of size: "<<j.first<<std::endl;
    // }

    /*=================================================================================
    ============================BEGIN CONSTRUCTING A SKETCH============================
    =================================================================================*/

    // initialize CM-sketch
    CM result(epsilon,delta);

    std::cout<<"Begin CM Sketch build"<<std::endl;
    omp_set_num_threads (numThreads);

    /* #pragma omp declare reduction(combineCM:CM: \
     omp_out=combineCM(omp_out,omp_in)) */

    #pragma omp parallel for shared(kmerLen,epsilon,delta,result) schedule(static,1)// reduction(combineCM:cm)
    for (unsigned int i=0;i<numThreads;i++) {
        CM cm(epsilon,delta);
        // cm.set_params(epsilon,delta);
        int curThread=omp_get_thread_num();
        // std::cout<<"THREAD #"<<curThread<<std::endl;
        for (auto fp : threadFPs[curThread]){
            // std::cout<<fp.second<<std::endl;
            std::string line, seq, c, cr;
            std::string curChrom="";
            std::ifstream genome(fp.second);
            if (genome.is_open()){
                while (std::getline(genome, line))
                {
                    if (line[0]=='>'){
                        curChrom=line.substr(1);
                        seq="";
                    }
                    else{
                        // now process kmers for each line and add to cmsketch
                        if (line.length()>=(unsigned int)kmerLen){ // shorter than kmerlength - typically last line of the fasta file
                            seq+=line.substr(0,kmerLen-1);
                        }
                        else{
                            seq+=line;
                        }
                        // std::cout<<line<<"\t"<<seq<<std::endl;
                        for (size_t k=0;k<seq.length()-kmerLen+1;k++){
                            // std::cout<<"forward: "<<seq.substr(k,kmerLen).c_str()<<std::endl;
                            c=seq.substr(k,kmerLen);
                            transform(c.begin(), c.end(), c.begin(), ::toupper);
                            cr=reverse(c);
                            if (cr.length()!=0){ //ambiguous code not detected
                                cm.update(c);
                                cm.update(cr);
                            }
                            // std::cout<<cm.estimate(seq.substr(k,kmerLen).c_str())<<std::endl;
                        }
                        if (line.length()>=(unsigned int)kmerLen){ // shorter than kmerlength - typically last line of the fasta file
                            for (size_t k=0;k<line.length()-kmerLen+1;k++){
                                // std::cout<<line.substr(k,kmerLen).c_str()<<std::endl;
                                c=line.substr(k,kmerLen);
                                transform(c.begin(), c.end(), c.begin(), ::toupper);
                                cr=reverse(c);
                                // std::cout<<c<<"\t"<<cr<<std::endl;
                                if (cr.length()!=0){ //ambiguous code not detected
                                    cm.update(c);
                                    cm.update(cr);
                                }
                            }
                            seq=line.substr(line.length()-kmerLen+1,kmerLen);
                        }
                        else{
                            seq=line;
                        }
                    }
                }
            }
            else{
                std::cout<<"file can not be opened"<<std::endl;
            }
            genome.close();
        }
        #pragma omp critical
        result.merge(&cm);
    }
    // SINGLE THREADED EVALUATION
    std::cout<<"Begin writing out results"<<std::endl;
    #pragma omp parallel for shared(kmerLen,result,args,max_count) schedule(static,1)// reduction(combineCM:cm)
    for (unsigned int i=0;i<numThreads;i++) {
        int curThread=omp_get_thread_num();
        int curCount=0;
        for (auto fp : threadFPs[curThread]){
            std::string line, seq, c;
            std::string curChrom="";
            std::ifstream genome(fp.second);
            FILE *outLog;
            std::string outFP=args.get_string(Opt::OUTPUT)+std::to_string(curThread);
            outLog = fopen(outFP.c_str(),"w+");

            if (genome.is_open()){
                long pos=0; // start of the current kmer
                while (std::getline(genome,line)){
                    if (line[0]=='>'){
                        curChrom=line.substr(1);
                        seq="";
                    }
                    else{
                        if (line.length()>=(unsigned int)kmerLen){ // shorter than kmerlength - typically last line of the fasta file
                            seq+=line.substr(0,kmerLen-1);
                        }
                        else{
                            seq+=line;
                        }
                        for (size_t k=0;k<seq.length()-kmerLen+1;k++){
                            c=seq.substr(k,kmerLen);
                            transform(c.begin(), c.end(), c.begin(), ::toupper);
                            curCount=result.estimate(c);
                            if (curCount<=max_count && curCount>0){
                                fprintf(outLog, "%s,%s,%ld,%d,%s\n",fp.second.c_str(),curChrom.c_str(),pos,curCount,c.c_str());
                            }
                            pos++;
                        }
                        if (line.length()>=(unsigned int)kmerLen){ // shorter than kmerlength - typically last line of the fasta file
                            for (size_t k=0;k<line.length()-kmerLen+1;k++){
                                c=line.substr(k,kmerLen);
                                transform(c.begin(), c.end(), c.begin(), ::toupper);
                                curCount=result.estimate(c);
                                if (curCount<=max_count && curCount>0){
                                    fprintf(outLog, "%s,%s,%ld,%d,%s\n",fp.second.c_str(),curChrom.c_str(),pos,curCount,c.c_str());
                                }
                                pos++;
                            }
                            seq=line.substr(line.length()-kmerLen+1,kmerLen);
                        }
                        else{
                            seq=line;
                        }
                    }
                }
            }
            fclose(outLog);
            genome.close();
        }
    }

    // now just merge all temporary files together

    return 0;
}