//[[Rcpp::depends(RcppArmadillo)]] 
#include<RcppArmadillo.h> 

using namespace Rcpp; 
using namespace arma; 

// [[Rcpp::export]] 
List JCDC(NumericMatrix Ar, NumericVector Similarity, int sim_len, NumericMatrix G0, double step_size, int K_fit, int iter_MAX, int grad_MAX, int Tabu_MAX, double MAX_weight, double lambda, int random_start=0, double alpha=1.0)
{
	int i, k, k1, z, count, converge;
	double step, old_reward, new_reward;
	
	mat A = as<mat>(Ar);
	int N = A.n_rows;

	cube phi(Similarity.begin(), N, N, sim_len);
	mat G_fit(G0.begin(), N, K_fit);
	mat theta(sim_len, K_fit), theta_new(sim_len, K_fit);
	
	mat Tabu_Method_2(N, K_fit);
	
	uword which_max;
	
	
	// 1. Check if G_fit is initialized correctly
	int degenerate=0;
	for(int check_k=0; check_k<K_fit; check_k++){
		if(sum(G_fit.col(check_k))<=1){degenerate=1;}
	}
	if(degenerate==1 || random_start==1){
		imat G_fit(N, K_fit);
		G_fit.zeros();
		
		for(i=0; i<N; i++){
			while(sum(G_fit.row(i))<0.5 || sum(G_fit.row(i))>2.5){
				for(k=0; k<K_fit; k++){
					 NumericVector temp1 = rbinom(1, 1, 1.0/K_fit);
					 G_fit(i, k) = temp1(0);
				}
			}
		}
	}
	
	// 2. Initialize parameters
	theta.zeros(); theta_new.zeros();
	
	// 3. Iterative estimation
	cube W(N, N, K_fit);
	mat W_scale(N, N), Ws_power(N, N);
	vec reward_vec(K_fit);
	vec gradient(sim_len);
	
	vec G_tmp0(N), G_tmp1(N);
	
	int Tabu_ind;
	//vec connectivity(K_fit);
	
	for(int iter_ind=0; iter_ind<iter_MAX; iter_ind++){
		
		// (I). Fix C, optimize over theta
		
		W_scale.zeros();
		for(k=0; k<K_fit; k++){
			Ws_power.zeros(); for(z=0; z<sim_len; z++){ Ws_power += theta(z, k)*phi.slice(z); }
			W.slice(k) = MAX_weight*ones(N, N) - trunc_exp( -Ws_power );
			if(sum(G_fit.col(k))>1.5){
				W_scale += (G_fit.col(k)*trans(G_fit.col(k)))/pow(sum(G_fit.col(k)), alpha) % W.slice(k);
			}
		}
		
		/////// calculate old_reward
		old_reward = accu( W_scale % A ) - lambda*accu(abs(theta));   //PENALIZED
		
		////// gradient
		for(k=0; k<K_fit; k++){
			gradient.zeros(); for(z=0; z<sim_len; z++){gradient(z) = accu( (G_fit.col(k)*trans(G_fit.col(k)))/pow(sum(G_fit.col(k)), alpha) % A % (MAX_weight*ones(N, N) - W.slice(k)) % phi.slice(z) ) - lambda*((theta(z, k)>0) - (theta(z, k)<0));}   //PENALIZED
			
			//// start gradient descent
			count = 0; converge = 0; step = step_size;
			while(count<grad_MAX && converge <5){
			
				theta_new.col(k) = theta.col(k) + step * gradient;
				for(z=0; z<sim_len; z++){
					if(theta_new(z, k)<0){
						theta_new(z, k) = 0;
					}
				}
				
				//// follow the same route to calculate new_loss
						W_scale.zeros();
						for(k1=0; k1<K_fit; k1++){
							Ws_power.zeros(); for(z=0; z<sim_len; z++){ Ws_power += theta_new(z, k1)*phi.slice(z); }
							W.slice(k1) = MAX_weight*ones(N, N) - trunc_exp( -Ws_power );
							if(sum(G_fit.col(k1))>1.5){
								W_scale += (G_fit.col(k1)*trans(G_fit.col(k1)))/pow(sum(G_fit.col(k1)), alpha) % W.slice(k1);
							}
						}
				
				/////// calculate new_loss
				new_reward = accu( W_scale % A ) - lambda*accu(abs(theta_new));   //PENALIZED
				
				if(new_reward>old_reward){
					old_reward = new_reward;
					theta = theta_new;
					count++;
					if( norm(gradient, 2) * step < 0.05 * norm(theta.col(k), 2) ){ converge++;	}	else{ converge = 0; }
					step *= 1.2;
					//////// update gradient; notice W is already updated
					gradient.zeros(); for(z=0; z<sim_len; z++){gradient(z) = accu( (G_fit.col(k)*trans(G_fit.col(k)))/pow(sum(G_fit.col(k)), alpha) % A % (MAX_weight*ones(N, N) - W.slice(k)) % phi.slice(z) ) - lambda*((theta(z, k)>0) - (theta(z, k)<0));}   //PENALIZED
				}
				else{
					count++;
					converge = 0;
					step *= 0.5;
				}
			
			}
			///
			///
			/// NOTICE THAT HERE WE DECIDED TO TRUNCATE ALL SOLVED THETA TO THEIR POSITIVE HALF - I.E. WE REQUIRE ALL DIMENSIONS OF THETA TO BE NONNEGATIVE.
			///
			///
			for(z=0; z<sim_len; z++){
				if(theta(z)<0){
					theta(z) = 0;
				}
			}
			
			W_scale.zeros();
			for(k1=0; k1<K_fit; k1++){
				Ws_power.zeros(); for(z=0; z<sim_len; z++){ Ws_power += theta(z, k1)*phi.slice(z); }
				W.slice(k1) = MAX_weight*ones(N, N) - trunc_exp( -Ws_power );
				if(sum(G_fit.col(k1))>1.5){
					W_scale += (G_fit.col(k1)*trans(G_fit.col(k1)))/pow(sum(G_fit.col(k1)), alpha) % W.slice(k1);
				}
			}
			
		}
		// end updating theta
		
		// (II). Fix theta, update C
		
		///// Notice that W has been updated to its best shape
		for(k=0; k<K_fit; k++){
			for(int ii=0; ii<N; ii++){
				for(int jj=0; jj<N; jj++){
					if(W(ii, jj, k)<0){
						W(ii, jj, k) = 0;
					}
				}
			}
		}// This truncates all negative weights in W into 0.	  

		for(Tabu_ind = 0; Tabu_ind<Tabu_MAX; Tabu_ind++){
			
			////// These two lines are used to accelerate - no need to recalculate Ebp once we fix theta.
			
			
			// Method 1:
			
			for(i=0; i<N; i++){
				reward_vec.zeros();
				for(k1=0; k1<K_fit; k1++){
				
					for(k=0; k<K_fit; k++){
						W_scale = W.slice(k);
						reward_vec(k) = sum(W_scale.col(i) % A.col(i) % G_fit.col(k)) / (pow(sum(G_fit.col(k))+0.01, alpha));
					}
					
				}
				
				G_fit.row(i) = trans(zeros(K_fit));
				reward_vec.max(which_max);
				G_fit(i, which_max) = 1;
				
			}
			
		
		}
	
	} // end big iteration
	
	return List::create(Named("G.fit")=G_fit, Named("theta")=theta, Named("step")=step, Named("gradient")=gradient, Named("W")=W);
}
