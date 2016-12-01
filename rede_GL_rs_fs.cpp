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
	//int buffersize;         	/* Number of steps with size h corresponding to time delay*/
	int layer;             		/* Neuron layer l2/3e=1, l2/3i=2 l4e=3, l4i=4, ..., l6e=7, l6i=8 */
	double v;                	/* Membrane potencial */
	double acum_ref;			/* Acumulator for the refractory time*/
	vector <double> w;                 	/* Synaptic weigth */	
	vector <int> delay;		/* Synapses delays */
	vector <double> syn_input;	/* List of synaptic increments to be added to the neuron considering time delay*/
	vector <int> post;        	/* Index of all post-synaptic neurons */
	vector <double> tdisp;    	/* Spiking times */
	double r,gamma,vth,delta,vs;	// Parameters GL-model
} ;

double * updateNeuronInput(int Ce, double rate, double J, double adjust_factor, double h, double tf){
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
    		input[(int)(values)] += (adjust_factor)*J;  // Excitatory input
    		randNum = gsl_rng_uniform_pos(r);
    		nextTime = -log(1-randNum) / aux_rate;
			values += nextTime;
	     	}
	}
	 
	input[0]=0;
	
	gsl_rng_free(r);
	
	return input;
}

/*
n[i].vth = vr + 0.5;	// Threshold potential. 0.5 is the multiplication of R*I, R = 10e-3 and I = 50.0
n[i].gamma = 0.037;	// gamma adjusted for RS Izhikevich neuron
n[i].delta = 0.0;	// delta adjusted for RS Izhikevich neuron
n[i].r = 1.0;		// r adjusted for RS Izhikevich neuron
n[i].vs = n[i].vth + pow((1-n[i].delta),(1/n[i].r))/n[i].gamma;
n[i].v = randun()*(n[i].vth-vr)+vr;
*/
void setParams(struct neuron *n, int type) {
	double vr = -65.0;
	// Ex	
	if(type == 0) {
		n->vth = vr + 0.5;
		n->gamma = 0.037;
		n->delta = 0.0;
		n->r = 1.0;
		n->vs = n->vth + pow((1-n->delta),(1/n->r))/n->gamma;
		n->v = randun()*(n->vth-vr)+vr;
	// In
	} else if(type == 1) {
		n->vth = vr + 1.92;
		n->gamma = 0.05;
		n->delta = 0.3;
		n->r = 1.0;
		n->vs = n->vth + pow((1-n->delta),(1/n->r))/n->gamma;
		n->v = randun()*(n->vth-vr)+vr;
	}
	
	 
}

void runsim(double W, double g, double seed) {

	/* W = Maximum synaptic weigth */
	/* g = Ratio between excitatory and inhibitory synaptic weights*/

	int net_size = 10000;				/* Number of neurons of the network*/
	double rescale_factor = net_size/80000.0;	/* Rescale factor of the network*/
	double h = 1.0;					/* Integration step */
	double tmax = 1000;          	/* Simulation time in ms */
	double de = 1.5;               	/* Delay for excitatory connections*/
	double di = 0.8;               	/* Delay for inhibitory connections */
	double t_ref = 2;				/* Refractory time*/
	double ** external_input = new double*[net_size];
		
	/* Neuron model parameters */
	double vr = -65.0;	// Resting Potencial
	double mi = 0.9;	// Decay constant
	double tau = 1.0;	// Time constant: IF model = 10; Izhi model = 1	
	
	double adjust_factor = tau/h; // Adjust factor for the EPSP: adjust_factor - Izhi model; tau/h - LIF model
	
	/* Using differents seeds to generate random numbers */
	srand (seed);
	/* Calculates the simulation time */
	clock_t startTime = clock();
	/* Number of simulation iterations */
	int n_it = round(tmax/h);

	/********************************************* CREATING NETWORK *******************************************************/

	cout << "Creating network: ";
	/* Creating all neurons */
	struct neuron* n = new struct neuron[net_size];

	/* Reading connections matrix */
	FILE *mat;
	mat = fopen("./ma10k.dat", "r");
	int pre = 0, pos = 0;
	while (fscanf(mat, "%d %d", &pre, &pos) != EOF) {
		n[pre].post.push_back(pos);  // Connect neurons by adding pre index to pres list of the postsynaptic neuron
	}
	fclose(mat);
	
	gsl_rng * r_aloc;
  	r_aloc = gsl_rng_alloc(gsl_rng_mt19937);
  	gsl_rng_set(r_aloc, chrono::high_resolution_clock::now().time_since_epoch().count());
  	
  	double delay_aux;
	int buffersize = ceil((de+0.75)/h)+1;

	/* Initiating variables for all neurons */
	for(int i = 0; i < net_size; i++) {
	  	/*
		n[i].vth = vr + 0.5;	// Threshold potential. 0.5 is the multiplication of R*I, R = 10e-3 and I = 50.0
		n[i].gamma = 0.037;	// gamma adjusted for RS Izhikevich neuron
		n[i].delta = 0.0;	// delta adjusted for RS Izhikevich neuron
		n[i].r = 1.0;		// r adjusted for RS Izhikevich neuron
		n[i].vs = n[i].vth + pow((1-n[i].delta),(1/n[i].r))/n[i].gamma;
		n[i].v = randun()*(n[i].vth-vr)+vr;
		*/
				
		buffersize = ceil((de+0.75)/h)+1;

		if(i>=0 && i < round(0.268* net_size)) {
			//external_input[i] = updateNeuronInput(2000, 8, W, 1, h, tmax);
			external_input[i] = updateNeuronInput(1600*rescale_factor, 8, W, adjust_factor, h, tmax);
			n[i].layer = 1; //L23e
			setParams(&n[i], 0);
			/* Initializing synaptic buffer with zeros*/
			for (unsigned int j = 0; j< n[i].post.size(); j++){
			  
				delay_aux=de + gsl_ran_gaussian(r_aloc,0.75);
				if(delay_aux<0) delay_aux=0;
				if(delay_aux>1.75) delay_aux=1.75;				
				
				n[i].delay.push_back(ceil(delay_aux/h));
				n[i].w.push_back((W + 0.1*gsl_ran_gaussian(r_aloc,W))*(adjust_factor));
			}
			
		} else if(i >= round(0.268* net_size) && i < round(0.344*net_size)) {
			//external_input[i] = updateNeuronInput(1850, 8, W, 1, h, tmax);
			external_input[i] = updateNeuronInput(1500*rescale_factor, 8, W, adjust_factor, h, tmax);
			n[i].layer = 2;	//L23i
			setParams(&n[i], 1);
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=di + gsl_ran_gaussian(r_aloc,0.4);
				if(delay_aux<0) delay_aux=0;
				if(delay_aux>1.75) delay_aux=1.75;				
				
				n[i].delay.push_back(ceil(delay_aux/h));
				n[i].w.push_back(-g*(W + 0.1*gsl_ran_gaussian(r_aloc,W))*(adjust_factor)); 
			}
			
		} else if(i >= round(0.344*net_size) && i < round(0.628*net_size))  {
			//external_input[i] = updateNeuronInput(2000, 8, W, 1, h, tmax);
			external_input[i] = updateNeuronInput(2100*rescale_factor, 8, W, adjust_factor, h, tmax);
			n[i].layer = 3;	//L4e
			setParams(&n[i], 0);
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=de + gsl_ran_gaussian(r_aloc,0.75);
				if(delay_aux<0) delay_aux=0;
				if(delay_aux>1.75) delay_aux=1.75;				
				
				n[i].delay.push_back(ceil(delay_aux/h));
				
				// Connections from L4e to L23e have synapse weight = 2*w
				if(n[n[i].post[j]].layer==1){	
					n[i].w.push_back(2*(W + 0.1*gsl_ran_gaussian(r_aloc,W))*(adjust_factor));
				}else{
					n[i].w.push_back((W + 0.1*gsl_ran_gaussian(r_aloc,W))*(adjust_factor));					
				}
			}
			
		} else if(i >= round(0.628*net_size) && i < round(0.700*net_size))  {
			//external_input[i] = updateNeuronInput(1850, 8, W, 1, h, tmax);
			external_input[i] = updateNeuronInput(1900*rescale_factor, 8, W, adjust_factor, h, tmax);
			n[i].layer = 4;	//L4i
			setParams(&n[i], 1);
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=di + gsl_ran_gaussian(r_aloc,0.4);
				if(delay_aux<0) delay_aux=0;
				if(delay_aux>1.75) delay_aux=1.75;				
				
				n[i].delay.push_back(ceil(delay_aux/h));
				n[i].w.push_back(-g*(W + 0.1*gsl_ran_gaussian(r_aloc,W))*(adjust_factor));
			}
			
		} else if(i >= round(0.700*net_size) && i < round(0.763*net_size))  {
			external_input[i] = updateNeuronInput(2000*rescale_factor, 8, W, adjust_factor, h, tmax);
			//external_input[i] = updateNeuronInput(2000, 8, W, 1, h, tmax);
			n[i].layer = 5;	//L5e
			setParams(&n[i], 0);
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=de + gsl_ran_gaussian(r_aloc,0.75);
				if(delay_aux<0) delay_aux=0;
				if(delay_aux>1.75) delay_aux=1.75;				
				
				n[i].delay.push_back(ceil(delay_aux/h));
				n[i].w.push_back((W + 0.1*gsl_ran_gaussian(r_aloc,W))*(adjust_factor));	 
			}
			
		} else if(i >= round(0.763*net_size) && i < round(0.777*net_size))  {
			external_input[i] = updateNeuronInput(1900*rescale_factor, 8, W, adjust_factor, h, tmax);
			//external_input[i] = updateNeuronInput(1850, 8, W, 1, h, tmax);
			n[i].layer = 6;	//L5i
			setParams(&n[i], 1);
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=di + gsl_ran_gaussian(r_aloc,0.4);
				if(delay_aux<0) delay_aux=0;
				if(delay_aux>1.75) delay_aux=1.75;				
				
				n[i].delay.push_back(ceil(delay_aux/h));
				n[i].w.push_back(-g*(W + 0.1*gsl_ran_gaussian(r_aloc,W))*(adjust_factor));	 
			}
			
		} else if(i >= round(0.777* net_size) && i < round(0.963*net_size)) {
			external_input[i] = updateNeuronInput(2900*rescale_factor, 8, W, adjust_factor, h, tmax);
			//external_input[i] = updateNeuronInput(2000, 8, W, 1, h, tmax);
			n[i].layer = 7;	//L6e
			setParams(&n[i], 0);
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=de + gsl_ran_gaussian(r_aloc,0.75);
				if(delay_aux<0) delay_aux=0;
				if(delay_aux>1.75) delay_aux=1.75;				
				
				n[i].delay.push_back(ceil(delay_aux/h));
				n[i].w.push_back((W + 0.1*gsl_ran_gaussian(r_aloc,W))*(adjust_factor));	 
			}
			
		} else if(i >= round(0.963*net_size) && i < round(1.0*net_size)) {
			external_input[i] = updateNeuronInput(2100*rescale_factor, 8, W, adjust_factor, h, tmax);
			//external_input[i] = updateNeuronInput(1850, 8, W, 1, h, tmax);
			n[i].layer = 8;	//L6i
			setParams(&n[i], 1);
			for (unsigned int j = 0; j< n[i].post.size(); j++){
				delay_aux=di + gsl_ran_gaussian(r_aloc,0.4);
				if(delay_aux<0) delay_aux=0;
				if(delay_aux>1.75) delay_aux=1.75;				
				
				n[i].delay.push_back(ceil(delay_aux/h));
				n[i].w.push_back(-g*(W + 0.1*gsl_ran_gaussian(r_aloc,W))*(adjust_factor));	 
			}
		}
		
		/* Initializing synaptic buffer with zeros*/
		for (int j = 0; j< buffersize; j++){
			n[i].syn_input.push_back(0);
		}
	}

	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " s.\n" << endl;
	startTime = clock();

	/********************************************************** SIMULATION **********************************************************/
	int ind_postlist = 0;		/* Index for the postsynaptic neuron*/
	double ext_input = vr;//-40;	/* External input*/
	int ind = 0;			/* Index for the buffer of synaptic inputs considering time delays*/
	double prob = 0.0;
	cout << "Simulation: " << endl;

	for(int i = 0; i < n_it; i++) {
		for(int j = 0; j < net_size; j++) {
			ind = (i)%(buffersize);

			//GL model
			if(n[j].acum_ref < t_ref ){
				n[j].acum_ref+=h;
			}else{
				n[j].v = mi * n[j].v + (1-mi) * vr + n[j].syn_input[ind] + external_input[j][i];
			}

			n[j].syn_input[ind] = 0;
			

			if(n[j].v >= n[j].vs) {
				n[j].tdisp.push_back(i);
				n[j].v = vr;
				n[j].acum_ref = 0;

				/* Adding synaptic inputs considering time delays between each pair of neurons connected*/
				for (unsigned int p = 0; p < n[j].post.size(); p++){
				    ind_postlist = n[j].post[p];
				    n[ind_postlist].syn_input[(ind + n[j].delay[p])%(buffersize)]+=n[j].w[p];
				}	
			}else if(n[j].v > n[j].vth && n[j].v < n[j].vs) {
				prob = pow(n[j].gamma*(n[j].v-n[j].vth), n[j].r) + n[j].delta;
				if (randun() < prob) {
				    n[j].tdisp.push_back(i);
				    n[j].v = vr;
    				    n[j].acum_ref = 0;

				    /* Adding synaptic inputs considering time delays between each pair of neurons connected*/
				    for (unsigned int p = 0; p < n[j].post.size(); p++){
					ind_postlist = n[j].post[p];
					n[ind_postlist].syn_input[(ind + n[j].delay[p])%(buffersize)]+=n[j].w[p];
				    }
				}
			}
		}
	}
	//cout << freq << endl;
	cout << "End of simulation" << endl << "Recording raster..." << endl;
	ostringstream nameFile;
	nameFile << "raster" /*<< net_size/1000 << "_w" << (int)round(W*1000) << "_g" << (int)(g*10)  << "_sd" << (int)seed */<<".dat";
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

	cout << "Time: "<< double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " s.\n\n" << endl;
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

	runsim(0.15,4,1);

	return 0;
}
