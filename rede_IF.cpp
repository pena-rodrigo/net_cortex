// to run: 

# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <vector>
# include <sys/time.h>
# include <ctime>
# include <time.h> 
//# include "mpi.h"
# include <sstream>
# include <string>
# include <fstream>
# include <unistd.h>
# include <algorithm>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
# include <chrono>

/* Generate a uniform randon number between 0 and 1 */
double randun() {
	return (double)rand() / (double)RAND_MAX;
}

using namespace std;

struct neuron {
	int buffersize;         	/* Number of steps with size h corresponding to time delay*/
	int layer;             		/* Neuron layer l2/3e=1, l2/3i=2 l4e=3, l4i=4, ..., l6e=7, l6i=8 */
	double v;                	/* Membrane potencial */
	double u;
	double acum_ref;			/* Acumulator for the refractory time*/
	vector <double> w;                 	/* Synaptic weigth */	
	vector <double> delay;		/* Synapses delays */
	vector <double> syn_input;	/* List of synaptic increments to be added to the neuron considering time delay*/
	vector <int> post;        	/* Index of all post-synaptic neurons */
	vector <double> tdisp;    	/* Spiking times */
	double a,b,c,d;	//parameters Izhikevich
} ;

double * updateNeuronInput(int Ce, double rate, double J, double tau, double h, double tf){
	double * input = new double[(int) floor(tf/h)+1];
	int N = ((int) floor(tf/h)+1);
	for (int i = 0; i < N; i++){
    	input[i] = 0;
    }
	//cout << "working" << endl;
  	gsl_rng *r;
  	double randNum;
  	double nextTime=0, values=0;

	if((r = gsl_rng_alloc(gsl_rng_mt19937)) == NULL) {
		printf("ERROR: Could not create random number generator\n");
	    exit(1);
	}
	
	gsl_rng_set(r, chrono::high_resolution_clock::now().time_since_epoch().count());
   
    // Poisson generator by method 1
 	double aux_rate = (rate*0.001*h);
	for (int n = 0; n < Ce; n++){
		values=0;	
    	while ((int)(values) < N){	
    		input[(int)(values)] += (1/h)*J*tau;  // Excitatory input
    		randNum = gsl_rng_uniform_pos(r);
    		nextTime = -log(1-randNum) / aux_rate;
			values += nextTime;
	     	}
	}
	 
	input[0]=0;
	
	gsl_rng_free(r);
	
	return input;
}

void runsim(double W, double g, double seed) {

	/* W = Maximum synaptic weigth */
	/* g = Ratio between excitatory and inhibitory synaptic weights*/

	int net_size = 4000;			/* Number of neurons of the network*/
	double h = 0.1;					/* Integration step */
	double tmax = 2000;          	/* Simulation time in ms */
	double de = 1.5;               	/* Delay for excitatory connections*/
	double di = 0.8;               	/* Delay for inhibitory connections */
	double t_ref = 2;				/* Refractory time*/
	double ** external_input = new double*[net_size];
	//for(int i = 0; i < net_size; i++){
//		 external_input[i] = new double[(int) floor(tmax/h)+1];
		// external_input[i] = updateNeuronInput(1000, 8, W, 1, h, tmax);
	//}
	
//	for (int j = 0; j < net_size; j++){
//   		for (int i = 0; i < (int) floor(tmax/h)+1; i++){	
//   			external_input[j][i] = 0;
//   		}
//   	}
//   	
//   	for (int j = 0; j < net_size; j++){
   		
   		
		
	/* Neuron model parameters */
	double vr = -65.0;           	/* Resting Potencial */
	double vth = 30;//-50.0;				/* Threshold potential*/
	double tau = 10.0;				/* Time constant for IF model */
	
	double a_FS=0.1,b_FS=0.2,c_FS=-65,d_FS=2,a_RS=0.02,b_RS=0.2,c_RS=-65,d_RS=2;
	
	/* Using differents seeds to generate random numbers */
	srand (seed);
	/* Calculates the simulation time */
	clock_t startTime = clock();
	/* Number of simulation iterations */
	int n_it = round(tmax/h);

	/********************************************************** CREATING NETWORK **********************************************************/

	cout << "Creating network: ";
	/* Creating all neurons */
	struct neuron* n = new struct neuron[net_size];

	/* Reading connections matrix */
	FILE *mat;
	mat = fopen("./ma4k.dat", "r");
	int pre = 0, pos = 0;
	while (fscanf(mat, "%d %d", &pre, &pos) != EOF) {
		n[pre].post.push_back(pos);  // Connect neurons by adding pre index to pres list of the postsynaptic neuron
	}
	fclose(mat);
	
	gsl_rng * r_aloc;
  	r_aloc = gsl_rng_alloc(gsl_rng_mt19937);
  	gsl_rng_set(r_aloc, chrono::high_resolution_clock::now().time_since_epoch().count());
  	
  	double delay_aux;
  	
	/* Initiating variables for all neurons */
	for(int i = 0; i < net_size; i++) {
		n[i].v = randun()*(vth-vr)+vr;
		n[i].u = n[i].b*n[i].v;
		n[i].buffersize = 3.5*round(de/h)+1;

		if(i>=0 && i < round(0.268* net_size)) {
		    external_input[i] = updateNeuronInput(2000, 8, W, 1, h, tmax);
			n[i].layer = 1; //L23e
			//n[i].ext_in = 1600*rf;
			/* Initializing synaptic buffer with zeros*/
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=de + gsl_ran_gaussian(r_aloc,0.75);
				if(delay_aux<0) delay_aux=0;
				
				n[i].delay.push_back(delay_aux);
				//cout << delay_aux << endl;
//				n[i].buffersize = round(n[i].delay.size())+1;
				n[i].w.push_back((W + 0.1*gsl_ran_gaussian(r_aloc,W))*(1/h));
				//cout << (W + 0.1*gsl_ran_gaussian(r_aloc,W)) << endl;
				
				n[i].a=a_RS; n[i].b=b_RS;	
				n[i].c=c_RS; n[i].d=d_RS;
			}
			
		} else if(i >= round(0.268* net_size) && i < round(0.344*net_size)) {
			external_input[i] = updateNeuronInput(1850, 8, W, 1, h, tmax);

			n[i].layer = 2;	//L23i
			//n[i].ext_in = 1500*rf;
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=di + gsl_ran_gaussian(r_aloc,0.4);
				if(delay_aux<0) delay_aux=0;
				
				n[i].delay.push_back(delay_aux);
//				n[i].buffersize = round(n[i].delay.size())+1;
				n[i].w.push_back(-g*(W + 0.1*gsl_ran_gaussian(r_aloc,W))*(1/h));
				
				n[i].a=a_FS; n[i].b=b_FS;	
				n[i].c=c_FS; n[i].d=d_FS;
			}
			
		} else if(i >= round(0.344*net_size) && i < round(0.628*net_size))  {
			external_input[i] = updateNeuronInput(2000, 8, W, 1, h, tmax);
			n[i].layer = 3;	//L4e
			//n[i].ext_in = 2100*rf;
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=de + gsl_ran_gaussian(r_aloc,0.75);
				if(delay_aux<0) delay_aux=0;
				
				n[i].delay.push_back(delay_aux);
//				n[i].buffersize = round(n[i].delay.size())+1;
				
				// Connections from L4e to L23e have synapse weight = 2*w
				if(n[n[i].post[j]].layer==1){	
					n[i].w.push_back(2*(W + 0.1*gsl_ran_gaussian(r_aloc,W))*(1/h));
				}else{
					n[i].w.push_back((W + 0.1*gsl_ran_gaussian(r_aloc,W))*(1/h));					
				}
								
				n[i].a=a_RS; n[i].b=b_RS;	
				n[i].c=c_RS; n[i].d=d_RS;
			}
			
		} else if(i >= round(0.628*net_size) && i < round(0.700*net_size))  {
			external_input[i] = updateNeuronInput(1850, 8, W, 1, h, tmax);
			n[i].layer = 4;	//L4i
			//n[i].ext_in = 1900*rf;
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=di + gsl_ran_gaussian(r_aloc,0.4);
				if(delay_aux<0) delay_aux=0;
				
				n[i].delay.push_back(delay_aux);
//				n[i].buffersize = round(n[i].delay.size())+1;
				n[i].w.push_back(-g*(W + 0.1*gsl_ran_gaussian(r_aloc,W))*(1/h));
				
				n[i].a=a_FS; n[i].b=b_FS;	
				n[i].c=c_FS; n[i].d=d_FS;
			}
			
		} else if(i >= round(0.700*net_size) && i < round(0.763*net_size))  {
		    external_input[i] = updateNeuronInput(2000, 8, W, 1, h, tmax);
			n[i].layer = 5;	//L5e
			//n[i].ext_in = 2000*rf;
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=de + gsl_ran_gaussian(r_aloc,0.75);
				if(delay_aux<0) delay_aux=0;
				
				n[i].delay.push_back(delay_aux);
//				n[i].buffersize = round(n[i].delay.size())+1;
				n[i].w.push_back((W + 0.1*gsl_ran_gaussian(r_aloc,W))*(1/h));
								
				n[i].a=a_RS; n[i].b=b_RS;	
				n[i].c=c_RS; n[i].d=d_RS;
			}
			
		} else if(i >= round(0.763*net_size) && i < round(0.777*net_size))  {
			external_input[i] = updateNeuronInput(1850, 8, W, 1, h, tmax);
			n[i].layer = 6;	//L5i
			//n[i].ext_in = 1900*rf;;
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=di + gsl_ran_gaussian(r_aloc,0.4);
				if(delay_aux<0) delay_aux=0;
				
				n[i].delay.push_back(delay_aux);
//				n[i].buffersize = round(n[i].delay.size())+1;
				n[i].w.push_back(-g*(W + 0.1*gsl_ran_gaussian(r_aloc,W))*(1/h));
				
				n[i].a=a_FS; n[i].b=b_FS;	
				n[i].c=c_FS; n[i].d=d_FS;
			}
			
		} else if(i >= round(0.777* net_size) && i < round(0.963*net_size)) {
			external_input[i] = updateNeuronInput(2000, 8, W, 1, h, tmax);
			n[i].layer = 7;	//L6e
			//n[i].ext_in = 2900*rf;
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=de + gsl_ran_gaussian(r_aloc,0.75);
				if(delay_aux<0) delay_aux=0;
				
				n[i].delay.push_back(delay_aux);
//				n[i].buffersize = round(n[i].delay.size())+1;
				n[i].w.push_back((W + 0.1*gsl_ran_gaussian(r_aloc,W))*(1/h));
								
				n[i].a=a_RS; n[i].b=b_RS;	
				n[i].c=c_RS; n[i].d=d_RS;
			}
			
		} else if(i >= round(0.963*net_size) && i < round(1.0*net_size)) {
			external_input[i] = updateNeuronInput(1850, 8, W, 1, h, tmax);
			n[i].layer = 8;	//L6i
			//n[i].ext_in = 2100*rf;
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=di + gsl_ran_gaussian(r_aloc,0.4);
				if(delay_aux<0) delay_aux=0;
				
				n[i].delay.push_back(delay_aux);
//				n[i].buffersize = round(n[i].delay.size())+1;
				n[i].w.push_back(-g*(W + 0.1*gsl_ran_gaussian(r_aloc,W))*(1/h));
								
				n[i].a=a_FS; n[i].b=b_FS;	
				n[i].c=c_FS; n[i].d=d_FS;
			}
		}
		
		/* Initializing synaptic buffer with zeros*/
		for (int j = 0; j< n[i].buffersize; j++){
			n[i].syn_input.push_back(0);
		}
	}

	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " s.\n" << endl;
	startTime = clock();

	/********************************************************** SIMULATION **********************************************************/
	int ind_postlist = 0;		/* Index for the postsynaptic neuron*/
	double ext_input = 0;//-40;	/* External input*/
	int ind = 0;			/* Index for the buffer of synaptic inputs considering time delays*/
	cout << "Simulation: " << endl;

	for(int i = 0; i < n_it; i++) {
		for(int j = 0; j < net_size; j++) {

			ind = (i)%(n[j].buffersize);
			
//			if (n[j].layer == 3){
		//		ext_input = -45;
//			}else{
//				ext_input = vr;
//			}

			if (i> 2000/h){
				ext_input = 0;
			}

			//IF model solved by Euler method
			if(n[j].acum_ref < t_ref ){
				n[j].acum_ref+=h;
			}else{
				//n[j].v += (h/tau)*(-n[j].v + ext_input + n[j].syn_input[ind]);
				n[j].v += h*(0.04*n[j].v*n[j].v + 5*n[j].v + 140 - n[j].u + ext_input + n[j].syn_input[ind] + external_input[j][i]);
				n[j].u += h*(n[j].a*(n[j].b*n[j].v - n[j].u));
			}

			n[j].syn_input[ind] = 0;

			if(n[j].v > vth){
				n[j].v = n[j].c;
				n[j].u += n[j].d;
				n[j].tdisp.push_back(double(i)*h);
				n[j].acum_ref=0;

				/* Adding synaptic inputs considering time delays between each pair of neurons connected*/
				for (unsigned int p = 0; p < n[j].post.size(); p++){
					ind_postlist = n[j].post[p];
					n[ind_postlist].syn_input[(ind + int(round((n[j].delay[p])/h)) - 1)%(n[j].buffersize)]+=n[j].w[p];
				}				
			}

		}
	}
	//cout << freq << endl;
	ostringstream nameFile;
	nameFile << "raster" << net_size/1000 << "_w" << (int)round(W*1000) << "_g" << (int)(g*10)  << "_sd" << (int)seed <<".dat";
	string str = nameFile.str();

	//chdir("data");
	FILE * data = fopen(str.c_str(), "w+");

	for(int i = 0; i < net_size; i++) {
		for(unsigned int j = 0; j < n[i].tdisp.size(); j++) {
			//if(n[i].tdisp[j]>=1000)
			fprintf(data, "%d\t%.1f\n", i, n[i].tdisp[j]);
		}
	}
	fclose(data);
	
	for(int i = 0; i < net_size; i++)
       	delete[] external_input[i];
    delete[] external_input;
    external_input = nullptr;
    
    gsl_rng_free(r_aloc);

	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " s.\n\n" << endl;
	return;
}

int main(int argc, char *argv[]){

	//compile with: mpic++ *.cpp
	//run with: mpirun -np NUMBEROFNODES a.out
	/*
	 *    int numTasks, rank, nsim = 30, ind = 0;
	 *    double params[nsim][4];
	 *
	 *    for (int i = 0; i < 10; i++){
	 *		for (int j = 0; j < 1; j++){
	 *			params[ind][0] = 0.019 + 0.0*(double(i));
	 *                        params[ind][1] = 4.0 + 0.0*(double(j));
	 *			params[ind][2] = i;
	 *			params[ind][3] = 0.019 + 0.00*(double(i));
	 *			ind ++;
}
}
	 */
	/*
	 *    params[0][0] = 0.15/8;
	 *    params[0][1] = 4.0;
	 *    params[0][2] = 3;
	 */

	/* MPI::Init(argc, argv);
	 *    numTasks = MPI::COMM_WORLD.Get_size();
	 *    rank = MPI::COMM_WORLD.Get_rank();
	 *
	 *
	 *	//for (int i = 0; i < nsim/numTasks + 1; i++){
	 *	//	ind = :q
	 *	//rank+numTasks*i;
	 *	//	if (ind < nsim){
	 *			cout << params[rank][0] << "\t" << params[rank][1] << "\t" << params[rank][2] << endl;
	 *			runsim(params[rank][0],params[rank][3],params[rank][1],params[rank][2]);
	 *	//	}
	 *	//}
	 *	MPI::Finalize();*/

	runsim(0.35,4,1);

	return 0;
}
