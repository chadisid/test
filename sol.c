#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef int32_t sox_int32_t;
typedef sox_int32_t sox_sample_t;
#define PI	3.14159265359
#define LSX_USE_VAR(x)  ((void)(x=0)) /* During static analysis, initialize unused variables to 0. */
#define SOX_INT_MIN(bits) (1 <<((bits)-1))
#define SOX_INT_MAX(bits) (((unsigned)-1)>>(33-(bits)))
#define SOX_SAMPLE_MAX (sox_sample_t)SOX_INT_MAX(32)
#define SOX_SAMPLE_TO_FLOAT_32BIT(d,clips) ((d)*(1.0 / (SOX_SAMPLE_MAX + 1.0)))
#define SOX_SAMPLE_MIN (sox_sample_t)SOX_INT_MIN(32)
#define SOX_SAMPLE_LOCALS sox_sample_t sox_macro_temp_sample; \
  double sox_macro_temp_double 
#define SOX_FLOAT_32BIT_TO_SAMPLE(d,clips) SOX_FLOAT_64BIT_TO_SAMPLE(d, clips)
#define SOX_FLOAT_64BIT_TO_SAMPLE(d, clips)                     \
  (sox_sample_t)(                                               \
    LSX_USE_VAR(sox_macro_temp_sample),                         \
    sox_macro_temp_double = (d) * (SOX_SAMPLE_MAX + 1.0),       \
    sox_macro_temp_double < 0 ?                                 \
      sox_macro_temp_double <= SOX_SAMPLE_MIN - 0.5 ?           \
        ++(clips), SOX_SAMPLE_MIN :                             \
        sox_macro_temp_double - 0.5 :                           \
      sox_macro_temp_double >= SOX_SAMPLE_MAX + 0.5 ?           \
        sox_macro_temp_double > SOX_SAMPLE_MAX + 1.0 ?          \
          ++(clips), SOX_SAMPLE_MAX :                           \
          SOX_SAMPLE_MAX :                                      \
        sox_macro_temp_double + 0.5                             \
  )
size_t clips_t = 0;
size_t clips_two = 0;
enum filter_type {
    highpass,
    lowpass,
};
typedef struct filter_context {
    double gain;
    double frequency;
    double mix;
    double i1, i2;
    double o1, o2;
    double a0, a1, a2;
    double b0, b1, b2;	
} filter_context;
void filter (float *input, float *output, int len,filter_context *filter)                                     
{                
    float *ibuf = input;                                                 
    float *obuf = output;                                                                                                         
    double wet = 1;                                      
    double dry = 1. - wet;                                                    
    double out;                                                               
    int i;                                                                    
    double a1 = -filter->a1;                                                                 
    double a2 = -filter->a2;                                                                 
    for (i = 0; i+1 < len; i++) {                                             
        filter->o2 = filter->i2 * filter->b2 + filter->i1 * filter->b1 + ibuf[i] * filter->b0 + filter->o2 * a2 + filter->o1 * a1;            
        filter->i2 = ibuf[i];                                                         
        out = filter->o2 * wet + filter->i2 * dry;                                                                                                        
        obuf[i] = out;                                                                                                                       
        i++;                                                                  
        filter->o1 = filter->i1 * filter->b2 + filter->i2 * filter->b1 + ibuf[i] * filter->b0 + filter->o1 * a2 + filter->o2 * a1;            
        filter->i1 = ibuf[i];                                                         
        out = filter->o1 * wet + filter->i1 * dry;                                                                                                        
        obuf[i] = out;                                                    
                                                                            
    }                                                                                      
        double o0 = ibuf[i] * filter->b0 + filter->i1 * filter->b1 + filter->i2 * filter->b2 + filter->o1 * a1 + filter->o2 * a2;     
        filter->i2 = filter->i1;                                                              
        filter->i1 = ibuf[i];                                                         
        filter->o2 = filter->o1;                                                              
        filter->o1 = o0;                                                              
        out = o0 * wet + filter->i1 * dry;                                                                                                        
        obuf[i] = out;                                                                                                                       
    }                                                                                                                                     
}

filter_context * init_filter(double frequency, int sample_rate, enum filter_type filter_name) {
   double w0 = 2 * PI * frequency / sample_rate;
   filter_context *filter = (filter_context *)malloc(sizeof(filter_context));
   memset(filter, 0, sizeof(filter_context));
   switch (filter_name) {
	case lowpass:	
	    filter->a0 = 1;
	    filter->a1 = -exp(-w0);
	    filter->a2 = 0;
	    filter->b0 = 1 + filter->a1;
	    filter->b1 = 0;
	    filter->b2 = 0;
	    break;
	case highpass:
	    filter->a0 = 1;
            filter->a1 = -exp(-w0);
            filter->a2 = 0;
            filter->b0 = (1 - filter->a1) / 2;
            filter->b1 = -filter->b0;
            filter->b2 = 0;   
	    break;
	default:
            free(filter);
            return NULL;
            break;		   
    }	
    return filter;	
}
int main(int argc, char **argv)
{
    const char *outfilename, *filename;
    FILE *f, *outfile;
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input file> <output file>\n", argv[0]);
        exit(1);
    }
    filename    = argv[1];
    outfilename = argv[2];
    f = fopen(filename, "rb");
    if (!f) {
        fprintf(stderr, "Could not open %s\n", filename);
        exit(1);
    }
    outfile = fopen(outfilename, "wb");
    if (!outfile) {
    	fprintf(stderr, "Could not open %s\n", outfilename);
        fclose(f);
        exit(1);
    }
    int64_t data_size, current_offset;
    fseek(f, 0, SEEK_END);
    data_size = ftell(f);
    fseek(f, 0, SEEK_SET);
    int num_samples = 1024;
    float *inbuf = (float *)malloc(sizeof(float)*num_samples);	
    float *out_lowpass = (float *)malloc(sizeof(float)*num_samples);
    float *out_highpass = (float*) malloc(sizeof(float)*num_samples);	
    int bytes; 
    double frequency_lowpass = 500;
    double frequency_highpass = 200;	
    int sample_rate = 48000;
    filter_context *filter_lowpass = init_filter(frequency_lowpass,sample_rate,lowpass); 
	
    filter_context *filter_highpass = init_filter(frequency_highpass,sample_rate,highpass); 	
    while (current_offset < data_size) {
        bytes = fread(inbuf, sizeof(float), num_samples, f);
        current_offset += sizeof(float)*bytes;
        int i = 0;
        for(i = 0; i<bytes; i++) {
	SOX_SAMPLE_LOCALS;
		sox_sample_t ty = SOX_FLOAT_32BIT_TO_SAMPLE(inbuf[i],clips_t);
		float td = SOX_SAMPLE_TO_FLOAT_32BIT(ty,clips_two);
		inbuf[i] = td;	 
        }
	filter (inbuf, out_highpass, bytes, filter_highpass);    
	filter (out_highpass, out_lowpass, bytes, filter_lowpass);
        for(i = 0; i<bytes; i++) {
	  out_lowpass[i] = 0.5*out_lowpass[i];
	}	
	fwrite(out_lowpass, sizeof(float), bytes, outfile);
            
     }
	printf("current offset %lld data_size %lld clips %i clips_two %i \n",current_offset, data_size, clips_t,clips_two);
	
	
    fflush(outfile);
    fclose(outfile);
    return 0;
}
