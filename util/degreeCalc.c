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
    std::ifstream file;
    file.open(argv[1]);
    std::string line;
    int lineNum = 0;
    int numDiag = 0;
    int* rows;
    int* cols;
    int matSize;
    while (std::getline(file, line))
    {   
        if(line.find('%') == std::string::npos)
        {
            int a, b;
            double c;
            std::stringstream lineStream(line);
            lineStream >> a >> b >> c;
            if(lineNum == 0)
            {
                matSize = a;
                rows = new int[a];
                cols = new int[b];
                for(int i = 0; i < a; i++)
                {
                    rows[i] = 0;
                    cols[i] = 0;
                }
            }
            else
            {
                rows[a-1]++;
                cols[b-1]++;
            }
            lineNum++;
        }
    }

    int min = INT32_MAX;
    int max = INT32_MIN;
    for(int i = 0; i < matSize; i++)
    {
        if(rows[i] < min)
        {
            min = rows[i];
        }
        if(cols[i] < min)
        {
            min = cols[i];
        }
        if(rows[i] > max)
        {
            max = rows[i];
        }
        if(cols[i] > max)
        {
            max = cols[i];
        }
    }
    file.close();
    delete [] rows;
    delete [] cols;   
    std::cout << min << " " << max << " " << std::endl;
}