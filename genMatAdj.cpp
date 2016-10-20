# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <sys/time.h>
# include <ctime>
# include <time.h> 

# define NN 4000

double randun() {
	return (double)rand() / (double)RAND_MAX;
}

using namespace std;

int main() {

	//system("clear");
	srand (time(NULL)*10000000);
	clock_t startTime = clock();

	int nsyn = 0,
        r = 0, 
	c = 0;
	
	if(1){ // turn on to Potdjans and Diesmann network
	
	
		double table[8][8] = 
		{ 
			{ 0.101, 0.169,	0.044,	0.082,	0.032,	0.0,	0.008,	0.0   },
			{ 0.135, 0.137,	0.032,	0.052,	0.075,	0.0,    0.004,	0.0   },
			{ 0.008, 0.006,	0.050,	0.135,	0.007,	0.0003,	0.045,	0.0   },
			{ 0.069, 0.003,	0.079,	0.160,	0.003,	0.0,	0.106,	0.0   },
			{ 0.100, 0.062,	0.051,	0.006,	0.083,	0.373,	0.020,	0.0   },
			{ 0.055, 0.027,	0.026,	0.002,	0.060,	0.316,	0.009,	0.0   },
			{ 0.016, 0.007,	0.021,	0.017,	0.057,	0.020,	0.040,	0.225 },
			{ 0.036, 0.001,	0.003,	0.001,	0.028,	0.008,	0.066,	0.144 },
		};
	
		cout << "Gerando matriz de adjacência... ";
		FILE *raster;
		raster = fopen("ma4k.dat","w");			
	    
	    for(int i = 0; i < NN; i++) {
	    	if(i >= 0 && i < round(268.0/1000 * NN)) {
	            c = 0;
		    }
	        if(i >= round(268.0/1000 * NN) && i < round(344.0/1000 * NN)) {
	            c = 1;
	        }
	        if(i >= round(344.0/1000 *NN) && i < round(628.0/1000 * NN)) {
	            c = 2;
	        }
	        if(i >= round(628.0/1000 * NN) && i < round(700.0/1000 * NN)){
	            c = 3;
	        }
	        if(i >= round(700.0/1000 * NN) && i < round(763.0/1000 * NN)) {
	            c = 4;
	        }
	        if(i >= round(763.0/1000 * NN) && i < round(777.0/1000 * NN)) {
	            c = 5;
	        }
	        if(i >= round(777.0/1000 * NN) && i < round(963.0/1000 * NN)) {
	            c = 6;
	        }
	        if(i >= round(963.0/1000 * NN) && i < round(1000.0 / 1000 * NN)) {
	            c = 7;
	        }    
	    	for(int j = 0; j < NN; j++) {
		   if(j >= 0 && j < round(268.0/1000 * NN)) {
	                if(randun() < table[0][c]) {
	                	fprintf(raster,"%d\t%d\n", i, j);
	                    nsyn++;
	                } 
	            }
	            if(j >= round(268.0/1000 * NN) && j < round(344.0/1000 * NN)) {
	                if(randun() < table[1][c]) {
	                	fprintf(raster,"%d\t%d\n", i, j);
	                    nsyn++;
	                }  
	            }
	            if(j >= round(344.0/1000 *NN) && j < round(628.0/1000 * NN)) {
	                if(randun() < table[2][c]) {
	                	fprintf(raster,"%d\t%d\n", i, j);
	                    nsyn++;
	                } 
	            }
	            if(j >= round(628.0/1000 * NN) && j < round(700.0/1000 * NN)){
	                if(randun() < table[3][c]) {
	                	fprintf(raster,"%d\t%d\n", i, j);
	                    nsyn++;
	                }  
	            }
	            if(j >= round(700.0/1000 * NN) && j < round(763.0/1000 * NN)) {
	                if(randun() < table[4][c]) {
	                	fprintf(raster,"%d\t%d\n", i, j);
	                    nsyn++;
	                } 
	            }
	            if(j >= round(763.0/1000 * NN) && j < round(777.0/1000 * NN)) {
	                if(randun() < table[5][c]) {
	                	fprintf(raster,"%d\t%d\n", i, j);
	                    nsyn++;
	                } 
	            }
	            if(j >= round(777.0/1000 * NN) && j < round(963.0/1000 * NN)) {
	                if(randun() < table[6][c]) {
	                	fprintf(raster,"%d\t%d\n", i, j);
	                    nsyn++;
	                } 
	            }
	            if(j >= round(963.0/1000 * NN) && j < round(1000.0 / 1000 * NN)) {
	                if(randun() < table[7][c]) {
	                	fprintf(raster,"%d\t%d\n", i, j);
	                    nsyn++;
	                } 
	            }   
	    	}
	    }
	
		cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " s.\n\n" << endl;
	        cout << "NN : " << NN << "  Nsyns : " << nsyn << "\n";

	}else{ // here homogeneous random network
	
//cout << "Gerando matriz de adjacência... ";
		FILE *raster;
		raster = fopen("random.dat","w");	
		int Ce=1000, sort; //excitatory and inhibitory inputs per neuron
		float gama=0.25;
		int Ni = (int) round(NN*gama);
		int Ci = (int) round(Ce*gama);
		
		
		for (int i = 0; i < (NN+Ni); i++) {  //exc
			for (int j = 0; j < Ce; j++) {
				sort = round(randun()*((NN)-1));
				while (sort==i){	
				      sort = round(randun()*((NN)-1));
				}	
				fprintf(raster,"%d\t%d\n", sort, i);
				nsyn++;
			}		
		}
			
		for (int i = 0; i < (NN+Ni); i++) {    //inh
			for (int j = 0; j < Ci; j++) {
				sort = round(randun()*Ni+NN);				
				while (sort==i){	
					sort = round(randun()*Ni+NN);
				}
				fprintf(raster,"%d\t%d\n", sort, i);
				nsyn++;
			}
		}
			
			
		cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " s.\n\n" << endl;
	        cout << "Ne + Ni : " << (NN+Ni) << "  Nsyns : " << nsyn << "\n";			
			
		fclose(raster);
	             	
			
	//	cout <<	round(randun()*NN) << endl;
			
		
	
	}


	return 0;	
}
