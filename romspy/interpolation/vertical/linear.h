extern void create_weights(float *__restrict__ depth_in_array, unsigned long long int depth_in_len,
                    float *__restrict__ depth_at_each_point_array,
                    unsigned long long int horizontal_slice_len, unsigned long long int depth_out_len,
                    float *__restrict__ weight_array);

extern void apply_weights(float *__restrict__ weight_arr, float *__restrict__ from_arr, float *__restrict__ to_arr,
                         unsigned long long int total_len);