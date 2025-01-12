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
    char data[100];
    std::ifstream source;
    source.open(argv[1]);
    std::string line;
    int lineNum = 0;
    int numDiag = 0;
    int* verticeMap = NULL;
    int a, b;
    double c;
    while (std::getline(source, line))
    {   
        if(line.find('%') == std::string::npos)
        {
            std::stringstream lineStream(line);
            lineStream >> a >> b >> c;
            if(lineNum == 0)
            {
                verticeMap = new int[a];
                break;
            }
        }
    }
    int counter = a;
    for(int i = 0; i < a; i++)
    {
        verticeMap[i] = counter;
        counter--;
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
                outfile << a << " " << b << " " << nnz << std::endl;    
            }
            else
            {
                lineStream >> a >> b >> c;
                outfile << verticeMap[a-1] << " " << verticeMap[b-1] << " " << c << std::endl; 
            }
            lineNum++;
        }
    }
    delete[] verticeMap;
    outfile.close();
    newSource.close();
    source.close();
}