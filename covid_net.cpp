#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <ctime> 
#include <cstdlib>
#include <stdio.h>  
using namespace std;

const float h = .5f; 	  
const int   Tsim = 10000000 / .5f; 
const int   Nexc = 100; 	 
const int   Ninh = 25;  	 
const int   Nneur = Nexc + Ninh;
const int   Ncon = Nneur*Nneur*0.1f; 
float Vms[Nneur][2]; 
float Ums[Nneur][2]; 
float Iex[Nneur]; 		
float Isyn[Nneur]; 	
int pre_conns[Ncon];   
int post_conns[Ncon]; 
float weights[Ncon]; 
float y[Ncon];   
float X[Ncon][2];
float Y[Ncon][2];	
float psc_excxpire_time = 4.0f; 
float extrasyn_excxpire_time = 100.0f; 
float astrocytic_relax_time = 80.0f; 
float X_thr = 5.6f; 
float astro_beta = 1.0f; 
float astro_gamma = 0.72f; 
float minWeight = 20.0f; 
float maxWeight = 60.0f; 
float virus_gamma = 0.0f; 
float Iex_max = 40.0f; 
float a = 0.02f;
float b = 0.5f;
float c = -40.0f; 
float d = 100.0f;
float k = 0.5f;
float Vr = -60.0f;
float Vt = -45.0f;
float Vpeak = 35.0f;  
float V0 = -60.0f; 
float U0 = 0.0f;  
float Cm = 50.0f;  
const int width = int(100 /h);
float sumSp[width];
int s = 0;
int SsumSp = 0;
float wsP[2];
int f= 0;
const int Num_exp = 1000;
float freq[Num_exp];

void init_connections(){
	for (int con_idx = 0; con_idx < Ncon;){
		int pre = rand() % Nneur;
		int post = rand() % Nneur;
		pre_conns[con_idx] = pre;
		post_conns[con_idx] = post;
		weights[con_idx] = (rand() % ((int)(maxWeight - minWeight) * 10)) / 10.0f + minWeight;
		if (pre >= Nexc){
			weights[con_idx] = -weights[con_idx];
		}
		y[con_idx] = 0.0f;
		X[con_idx][0] = 0.0f;
		Y[con_idx][0] = 0.0f;
		X[con_idx][1] = 0.0f;
		Y[con_idx][1] = 0.0f;
		con_idx++;
	}
}

void init_neurons(){
	for (int neur_idx = 0; neur_idx < Nneur; neur_idx++){
		Iex[neur_idx] = ((float)rand() / RAND_MAX) * Iex_max;
		Isyn[neur_idx] = 0.0f;
		Vms[neur_idx][0] = V0;
		Ums[neur_idx][0] = U0;
		Vms[neur_idx][1] = 0;
		Ums[neur_idx][1] = 0;
		}
}


float izhik_Vm(int neuron){
	return (k*(Vms[neuron][0] - Vr)*(Vms[neuron][0] - Vt) - Ums[neuron][0] + Iex[neuron] + Isyn[neuron]) / Cm;
}

float izhik_Um(int neuron){
	return a*(b*(Vms[neuron][0] - Vr) - Ums[neuron][0]);
}

float astrocytic_glu(int con){
	return (- (1/ astrocytic_relax_time) * Y[con][0]  + (astro_beta * (1 - virus_gamma)) / (1 + exp(-X[con][0] + X_thr)));
}

void init_spike_window(){
		for (int iter = 0; iter < width; iter++){
			sumSp[iter] = 0;
		}
}

float spike_window(){
		if (s < width){
				for (int neur = 0; neur < Nneur; neur++){
					if (Vms[neur][0] >Vpeak){
						sumSp[s] +=1;
					}
				}
				s +=1;
			} else{
				for (int iter = 0; iter < width; iter++){
					SsumSp += sumSp[iter];
				}
				for (int iter = 1; iter < width; iter++){
					sumSp[iter - 1] = sumSp[iter];
					s = iter;
					sumSp[iter] = 0;
				}
				wsP[1] = SsumSp;

				if (wsP[0] >= 65 and wsP[1] < 65){
					f +=1;
				}

				wsP[0]= wsP[1];

				SsumSp=0;
			}
			return wsP[0];
}

int main(){

srand (1);

float sum_freq = 0.0f;

for (int m = 0; m < 100; m +=1){

sum_freq = 0.0f;

virus_gamma = float(m)/100;

	s = 0;
  SsumSp = 0;
  wsP[0] =0.0f;
	wsP[1] =0.0f;
  f = 0;

	init_connections();
	init_neurons();
	init_spike_window();

	float expire_coeff = exp(-h / psc_excxpire_time);
	float extrasynaptic_expire_coeff = exp(-h / extrasyn_excxpire_time);

ofstream res_file1;
ofstream res_file2;

	float T=0.0f;
	float T2=0.0f;

std::ostringstream fileNameStream1;
fileNameStream1 << "rastr_" << m << ".csv";
std::string fileName1 = fileNameStream1.str();
res_file1.open(fileName1.c_str());

std::ostringstream fileNameStream2;
fileNameStream2 << "oscill_" << m << ".csv";
std::string fileName2 = fileNameStream2.str();
res_file2.open(fileName2.c_str());

	for (int t = 1; t < Tsim; t++){
		float Vm_mean = 0.0f;
					
		for (int neur = 0; neur < Nneur; neur++){
			Vms[neur][1] = Vms[neur][0] + h*izhik_Vm(neur);
			Ums[neur][1] = Ums[neur][0] + h*izhik_Um(neur);
			Isyn[neur] = 0.0f;

			if (Vms[neur][0] >Vpeak){
				Vms[neur][1] = c;
				Ums[neur][1] = Ums[neur][0] + d;
				res_file1 << t*h << "; " << neur + 1 << "; " << endl;

			} 

			Vm_mean += Vms[neur][1];
		}

			Vm_mean /= Nneur;

			spike_window();

		for (int con = 0; con < Ncon; con++){
			y[con] = y[con] * expire_coeff;
			X[con][1] = X[con][0] * extrasynaptic_expire_coeff;
			Y[con][1] = Y[con][0] + h * astrocytic_glu(con);

			if (Vms[pre_conns[con]][0] > Vpeak){
				y[con] += 1.0f;
				X[con][1] += 1.0f;
			}

			if (weights[con] > 0){

			Isyn[post_conns[con]] += y[con] * weights[con] * (1 + astro_gamma * Y[con][1]);
			} else{
			Isyn[post_conns[con]] += y[con] * weights[con];

			}

			X[con][0] = X[con][1];
			Y[con][0] = Y[con][1];

		}



	//	if((t*h) > T){
		
				res_file2 << t*h << "; "  << Vm_mean <<"; " << endl;

		
	//			T = t*h + 0.1; //0.01
	//		}
		
		for (int neur = 0; neur < Nneur; neur++){
			Vms[neur][0] = Vms[neur][1];
			Ums[neur][0] = Ums[neur][1];
		}

	}

res_file1.close();
res_file2.close();	
	
}



	return 0;

}