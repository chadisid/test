#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef int32_t sox_int32_t;
typedef sox_int32_t sox_sample_t;
#define PI	3.14159265359
#define LSX_USE_VAR(x)  ((void)(x=0)) /* During static analysis, initialize unused variables to 0. */
/**
Client API:
Max value for sox_sample_t = 0x7FFFFFFF.
*//**
Client API:
Returns the smallest (negative) value storable in a twos-complement signed
integer with the specified number of bits, cast to an unsigned integer;
for example, SOX_INT_MIN(8) = 0x80, SOX_INT_MIN(16) = 0x8000, etc.
@param bits Size of value for which to calculate minimum.
@returns the smallest (negative) value storable in a twos-complement signed
integer with the specified number of bits, cast to an unsigned integer.
*/
#define SOX_INT_MIN(bits) (1 <<((bits)-1))

/**
Client API:
Returns the largest (positive) value storable in a twos-complement signed
integer with the specified number of bits, cast to an unsigned integer;
for example, SOX_INT_MAX(8) = 0x7F, SOX_INT_MAX(16) = 0x7FFF, etc.
@param bits Size of value for which to calculate maximum.
@returns the largest (positive) value storable in a twos-complement signed
integer with the specified number of bits, cast to an unsigned integer.
*/
#define SOX_INT_MAX(bits) (((unsigned)-1)>>(33-(bits)))


#define SOX_SAMPLE_MAX (sox_sample_t)SOX_INT_MAX(32)
#define SOX_SAMPLE_TO_FLOAT_32BIT(d,clips) ((d)*(1.0 / (SOX_SAMPLE_MAX + 1.0)))

/**
Client API:
Min value for sox_sample_t = 0x80000000.
*/
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
/* Helper struct to generate and parse WAVE headers */
typedef struct wav_header {
	char riff[4];
	uint32_t len;
	char wave[4];
	char fmt[4];
	uint32_t formatsize;
	uint16_t format;
	uint16_t channels;
	uint32_t samplerate;
	uint32_t avgbyterate;
	uint16_t samplebytes;
	uint16_t channelbits;
	char data[4];
	uint32_t blocksize;
} wav_header;
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

static void filter (float *input, float *output, int len,filter_context *filter)                                     
{    
    //len is number of samples	
    printf("Inside biquad filter \n");                   
    float *ibuf = input;                                                 
    float *obuf = output;                                                                                                         
    double wet = 1; //s->mix = 1;                                                      
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
    if (i < len) {                                                            
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
            //av_assert0(0);
            break;		   
    }	
    return filter;	
}
int main(int argc, char **argv)
{
    const char *outfilename, *filename;
    FILE *f, *outfile, *outraw;
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
        exit(1);
    }
    outraw = fopen("xx.raw", "wb");
    if (!outraw) {
        fprintf(stderr, "Could not open\n" );
        exit(1);
    }
    /* Write WAV header */
    wav_header header = {
	    {'R', 'I', 'F', 'F'},
	    0,
	    {'W', 'A', 'V', 'E'},
	    {'f', 'm', 't', ' '},
	    16,
	    0x0003,
	    1,
	    48000,
	    48000 * 2,
	    2,
	    16,
	    {'d', 'a', 't', 'a'},
	    0
	    };
    if(fwrite(&header, 1, sizeof(header), outfile) != sizeof(header)) {
	printf("Error writing WAV header...\n");
    }
    fflush(outfile);
    int64_t data_size, current_offset;
    fseek(f, 0, SEEK_END);
    data_size = ftell(f);
    printf("Size of input file %lld\n",data_size);
    fseek(f, 0, SEEK_SET);
    float max_val,min_val,max_val_two,min_val_two;
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
		if(inbuf[i] > max_val)  {
		    max_val = inbuf[i];
		}    
		if(inbuf[i] < min_val) {
		    min_val = inbuf[i];
		}
		//inbuf[i] = 0.05*inbuf[i];
		SOX_SAMPLE_LOCALS;
		sox_sample_t ty = SOX_FLOAT_32BIT_TO_SAMPLE(inbuf[i],clips_t);
		float td = SOX_SAMPLE_TO_FLOAT_32BIT(ty,clips_two);
		inbuf[i] = td;
		if(td > max_val_two)  {
		    max_val_two = td;
		}
		if(td < min_val_two) {
		    min_val_two = td;
		}
        }
	filter (inbuf, out_highpass, bytes, filter_highpass);    
	filter (out_highpass, out_lowpass, bytes, filter_lowpass);
        for(i = 0; i<bytes; i++) {
	  out_lowpass[i] = 0.5*out_lowpass[i];
	}	
	fwrite(out_lowpass, sizeof(float), bytes, outfile);
        fwrite(out_lowpass, sizeof(float), bytes, outraw);
        printf("Sample values %f current offset %lld data_size %lld max_val %f min_val %f clips %i clips_two %i min_val_two %f max_val_two %f\n",inbuf[0],current_offset, data_size, max_val,min_val,clips_t,clips_two,min_val_two,max_val_two);    
     }
     fseek(outfile, 0, SEEK_END);
     long int size = ftell(outfile);
     if(size >= 8) {
	size -= 8;
	fseek(outfile, 4, SEEK_SET);
	fwrite(&size, sizeof(uint32_t), 1, outfile);
	size += 8;
	fseek(outfile, 40, SEEK_SET);
	fwrite(&size, sizeof(uint32_t), 1, outfile);
	fflush(outfile);
	fseek(outfile, 0, SEEK_END);
    }
    fclose(outfile);
    fflush(outraw);
    fclose(outraw);
    return 0;
}
