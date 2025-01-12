#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>
#include <fstream>
#include <iostream>
#include <sstream>



int main(int argc, char *argv[]) 
{
    if(argc < 3)
    {
        std::cout << "Two file file arguments expected." << std::endl;
    }
    char data[100];
    std::ifstream source;
    source.open(argv[1]);
    std::string line;
    int lineNum = 0;
    int numDiag = 0;
    std::cout << "Mirroring edges for " << argv[1] << " and writing to " << argv[2] << std::endl;
    while (std::getline(source, line))
    {   
        if(line.find('%') == std::string::npos)
        {
            int a, b;
            double c;
            std::stringstream lineStream(line);
            lineStream >> a >> b >> c;
            if(lineNum != 0)
            {
                if(a == b)
                {
                    numDiag++;
                }
            }
            lineNum++;
        }
        
    }
    
    std::ifstream newSource;
    std::ofstream outfile (argv[2], std::ofstream::out);
    newSource.open(argv[1]);
    lineNum = 0;
    while (std::getline(newSource, line))
    {   
        if(line.find('%') == std::string::npos)
        {
            int a, b;
            double c;
            int nnz;
            std::stringstream lineStream(line);
            if(lineNum == 0)
            {
                lineStream >> a >> b >> nnz;
                outfile << a << " " << b << " " << (nnz * 2) - numDiag << " " << std::endl;    
            }
            else
            {
                lineStream >> a >> b >> c;
                outfile << a << " " << b << " " << c << " " << std::endl; 
                if(a != b)
                {
                    outfile << b << " " << a << " " << c << " " << std::endl; 
                }
            }
            lineNum++;
        }
    }
    
    source.close();
    newSource.close();
    outfile.close();
}