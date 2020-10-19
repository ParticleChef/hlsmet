/*
  HLS implementation of MET calculation from PF objects
*/
#include "met.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#include <string>
#endif

// pt, phi are integers
void met_hw(hls::stream<PFInputWords> &inputs, hls::stream<word_t> &output){
//#pragma HLS dataflow
#pragma HLS pipeline ii=36
    //#pragma HLS ARRAY_PARTITION variable=inputs complete

	if(DEBUG) std::cout << " ==== HW Begin ====" << std::endl;
    pt_t data_pt[NPART]; phi_t data_phi[NPART];  //var_t is ap_int<16>

    PFInputWords input_array = inputs.read();
    
    for(int i = 0; i < NPART; i++){
        data_pt[i]  = input_array.data[i](63,48);
        data_phi[i] = input_array.data[i](47,32);
        if(DEBUG2) std::cout<<"    >> input pt phi: "<<data_pt[i]<<" "<<data_phi[i]<<std::endl;
    }
#pragma HLS ARRAY_PARTITION variable=data_pt complete
#pragma HLS ARRAY_PARTITION variable=data_phi complete

    // calc signed components first
    pxy_t sum_x = 0;
    pxy_t sum_y = 0;

 LOOP_COMPONENTS: for(int i=0; i<NPART;i++){
        // Get x, y components
        pxy_t met_x = 0;
        pxy_t met_y = 0;
        ProjX(data_pt[i], data_phi[i], met_x);
        ProjY(data_pt[i], data_phi[i], met_y);
        // Add to x, y sums
        sum_x -= met_x;
        sum_y -= met_y;
        if(DEBUG2){
            std::cout << "    met x,y = (" << -met_x << ", " << -met_y << ") \t";
            std::cout << " sum x,y = (" << sum_x << ", " << sum_y << ") \t";
            std::cout << " \n";
        }
    }

    pt2_t res_pt2;
    res_pt2 = sum_x*sum_x + sum_y*sum_y;

    pt_t res_pt; phi_t res_phi;
    PhiFromXY(sum_x,sum_y,res_pt2,res_phi);

    res_pt = hls::sqrt(res_pt2);

    if(DEBUG) std::cout<<"    hw out res_pt&phi : "<<res_pt<<", "<<res_phi<<std::endl;
    if(DEBUG) std::cout<<"    hw conv res_pt&phi : "<<res_pt/4<<", "<<res_phi/(4*(180/3.141))<<std::endl;
	if(DEBUG) std::cout<<std::endl;

    word_t outword = ( res_pt, res_phi, ap_uint<32>(0) );
    output.write( outword );
    return;


    return;
}

