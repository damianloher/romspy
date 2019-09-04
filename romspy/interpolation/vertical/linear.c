#include <omp.h>
#include <stdio.h>

/*
 * Functions to interpolate vertically with weights.
 * Arrays have the dimension order (z, y, x).
 * Weights are of length 3 times the target 3D array.
 * To change interpolation method, write new interpolation method and increase weight length both here and in roms.
 * Author: Nicolas Munnich
 * License: GNU GPL2+
 * */


void create_weights(float *__restrict__ depth_in_array, unsigned long long int depth_in_len,
                    float *__restrict__ depth_at_each_point_array,
                    unsigned long long int horizontal_slice_len, unsigned long long int depth_out_len,
                    float *__restrict__ weight_array) {
    /*
     * Calculates weights for a 1D interpolation on a 3D array.
     * depth_in_array : input 1D array of linear depths (positive!)
     * depth_in_len : length of depth_in_array
     * depth_at_each_point_array : 3D array describing target depths. Must be (depth, eta, xi)
     * horizontal_slice_len : length of xi axis times length of eta axis
     * depth_out_len : length of depth out axis
     * weight_array : 4D array where the first arrays are of the same size as depth_at_each_point_array
     *                and the last array is a pair of [lower index, upper index,  percent]
     * */
    if (depth_in_array[depth_in_len >> 1] > 0) {
        for (unsigned int i = 0; i < depth_in_len; i++) {
            depth_in_array[i] = -depth_in_array[i];
        }
    }
#pragma omp parallel
    {
        float target_depth_value;
        unsigned int depth_in_index = 0;
        float percent;
        unsigned int absolute_index;

#pragma omp for
        for (unsigned int depth = 0; depth < depth_out_len; depth++) {
            for (unsigned int index = 0; index < horizontal_slice_len; index++) {

                absolute_index = depth * horizontal_slice_len + index;
                target_depth_value = depth_at_each_point_array[absolute_index];

                /*
                 * [0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10]
                 * -4.5
                 *
                 * */
                //find corresponding indexes
                if (depth_in_array[depth_in_index] < target_depth_value) { // If h_in
                    while (depth_in_index > 0) {
                        if (depth_in_array[depth_in_index] >
                            target_depth_value) { // shift until depth_in_index is one above h then subtract one
                            break;
                        }
                        depth_in_index--;
                    }
                } else if (depth_in_array[depth_in_index + 1] >= target_depth_value) { // If depth_in_index is above h
                    while (depth_in_index + 1 < depth_in_len - 1) {
                        if (depth_in_array[depth_in_index + 1] < target_depth_value) {
                            break;
                        }
                        depth_in_index++;

                    }
                }
                // Calculate corresponding percent
                percent = 1 - (target_depth_value - depth_in_array[depth_in_index]) /
                              (depth_in_array[depth_in_index + 1] - depth_in_array[depth_in_index]);
                if (percent > 1) {
                    percent = 1;
                } else if (percent < 0) {
                    percent = 0;
                }


                // Insert both indexes to increase of the interpolation at cost to memory
                weight_array[(absolute_index * 3)] = (float) (depth_in_index * horizontal_slice_len + index);
                weight_array[(absolute_index * 3) + 1] = (float) ((depth_in_index + 1) * horizontal_slice_len + index);
                weight_array[(absolute_index * 3) + 2] = percent;


            }
        }
    }
}

void apply_weights(float *__restrict__ weight_arr, float *__restrict__ from_arr, float *__restrict__ to_arr,
                         unsigned long long int total_len) {
    /* 
     * Interpolates a 3D array. The first index is the dimension which is interpolated. Uses weights.
     * weight_arr: 4D array containing weights
     * from_arr: 3D array containing inputs
     * to_arr: 3D array where output is placed
     * total_len: number of elements in to_arr
     * */
#pragma omp parallel
    {
        unsigned long long int index1;
        unsigned long long int index2;
        float index_percent;

#pragma omp for
        for (unsigned long long int index = 0; index < total_len; index++) {
            index1 = (unsigned long long int) weight_arr[(index * 3)]; // lower index
            index2 = (unsigned long long int) weight_arr[(index * 3) + 1]; // upper index
            index_percent = weight_arr[(index * 3) + 2]; // percentage lower index
            //Interpolate
            to_arr[index] = (float) (from_arr[index1] * index_percent + from_arr[index2] * (1 - index_percent));
        }
    }
}

//void interpolate(float *__restrict__ depth_in_array, unsigned long long int depth_in_len,
//                    float *__restrict__ depth_at_each_point_array,
//                    unsigned long long int horizontal_slice_len, unsigned long long int depth_out_len,
//                    float *__restrict__ from_arr, float *__restrict__ to_arr)){}