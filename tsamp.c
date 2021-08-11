#include <string.h>

#include <stdio.h>

#include <stdlib.h>

#include <stdint.h>

#include <math.h>

#define PI 3.14159265359
#define SAMPLE_MAX (int32_t)(1 <<31)
#define SAMPLE_MIN (int32_t)(((unsigned)-1)>>1)

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
}
filter_context;
void filter(float * input, float * output, int len, filter_context * filter) {
  float * ibuf = input;
  float * obuf = output;
  double wet = 1;
  double dry = 1. - wet;
  double out;
  int i;
  double a1 = -filter -> a1;
  double a2 = -filter -> a2;
  for (i = 0; i + 1 < len; i++) {
    filter -> o2 = filter -> i2 * filter -> b2 + filter -> i1 * filter -> b1 + ibuf[i] * filter -> b0 + filter -> o2 * a2 + filter -> o1 * a1;
    filter -> i2 = ibuf[i];
    out = filter -> o2 * wet + filter -> i2 * dry;
    obuf[i] = out;
    i++;
    filter -> o1 = filter -> i1 * filter -> b2 + filter -> i2 * filter -> b1 + ibuf[i] * filter -> b0 + filter -> o1 * a2 + filter -> o2 * a1;
    filter -> i1 = ibuf[i];
    out = filter -> o1 * wet + filter -> i1 * dry;
    obuf[i] = out;

  }
  if (i < len) {
  double o0 = ibuf[i] * filter -> b0 + filter -> i1 * filter -> b1 + filter -> i2 * filter -> b2 + filter -> o1 * a1 + filter -> o2 * a2;
  filter -> i2 = filter -> i1;
  filter -> i1 = ibuf[i];
  filter -> o2 = filter -> o1;
  filter -> o1 = o0;
  out = o0 * wet + filter -> i1 * dry;
  obuf[i] = out;
}
}

filter_context * init_filter(double frequency, int sample_rate, enum filter_type filter_name) {
  double w0 = 2 * PI * frequency / sample_rate;
  filter_context * filter = (filter_context * ) malloc(sizeof(filter_context));
  memset(filter, 0, sizeof(filter_context));
  switch (filter_name) {
  case lowpass:
    filter -> a0 = 1;
    filter -> a1 = -exp(-w0);
    filter -> a2 = 0;
    filter -> b0 = 1 + filter -> a1;
    filter -> b1 = 0;
    filter -> b2 = 0;
    break;
  case highpass:
    filter -> a0 = 1;
    filter -> a1 = -exp(-w0);
    filter -> a2 = 0;
    filter -> b0 = (1 - filter -> a1) / 2;
    filter -> b1 = -filter -> b0;
    filter -> b2 = 0;
    break;
  default:
    free(filter);
    return NULL;
    break;
  }
  return filter;
}
int main(int argc, char ** argv) {
  const char * outfilename, * filename;
  FILE * f, * outfile;
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <input file> <output file>\n", argv[0]);
    exit(1);
  }
  filename = argv[1];
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
  float * inbuf = (float * ) malloc(sizeof(float) * num_samples);
  if (!inbuf) {
    fprintf(stderr, "Could not allocate memory\n");
    fclose(outfile);
    fclose(f);
    exit(1);
  }
  float * out_lowpass = (float * ) malloc(sizeof(float) * num_samples);
  if (!out_lowpass) {
    fprintf(stderr, "Could not allocate memory\n");
    free(inbuf);
    fclose(outfile);
    fclose(f);
    exit(1);
  }
  float * out_highpass = (float * ) malloc(sizeof(float) * num_samples);
  if (!out_highpass) {
    fprintf(stderr, "Could not allocate memory\n");
    free(out_lowpass);
    free(inbuf);
    fclose(outfile);
    fclose(f);
    exit(1);
  }
  int bytes;
  double frequency_lowpass = 500;
  double frequency_highpass = 200;
  int sample_rate = 48000;
  filter_context * filter_lowpass = init_filter(frequency_lowpass, sample_rate, lowpass);
  if (!filter_lowpass) {
    fprintf(stderr, "Could not initialize low pass filter\n");
    free(out_highpass);
    free(out_lowpass);
    free(inbuf);
    fclose(outfile);
    fclose(f);
    exit(1);
  }
  filter_context * filter_highpass = init_filter(frequency_highpass, sample_rate, highpass);
  if (!filter_highpass) {
    fprintf(stderr, "Could not initialize high pass filter\n");
    free(filter_lowpass);
    free(out_highpass);
    free(out_lowpass);
    free(inbuf);
    fclose(outfile);
    fclose(f);
    exit(1);
  }
  size_t clips_on = 0;
  while (current_offset < data_size) {
    bytes = fread(inbuf, sizeof(float), num_samples, f);
    current_offset += sizeof(float) * bytes;
    int i;
    for (i = 0; i < bytes; i++) {
      int32_t temp_sample_test;
      double sox_macro_temp_double_temp = (inbuf[i]) * (SAMPLE_MAX + 1.0);
      if(sox_macro_temp_double_temp < 0) {
         if(sox_macro_temp_double_temp <= SAMPLE_MIN - 0.5){
          ++(clips_on);
          temp_sample_test= SAMPLE_MIN ;
         } else {
          temp_sample_test = sox_macro_temp_double_temp - 0.5;
         }
      } else {
         if(sox_macro_temp_double_temp >= SAMPLE_MAX + 0.5 ){
            if( sox_macro_temp_double_temp > SAMPLE_MAX + 1.0 ){
              ++(clips_on);
              temp_sample_test = SAMPLE_MAX ;
            } else {
              temp_sample_test = SAMPLE_MAX ;
            }
         } else {
           temp_sample_test = sox_macro_temp_double_temp + 0.5;
         }
      }
      float temp_sample_float_t = (temp_sample_test)*(1.0 / (SAMPLE_MAX + 1.0));
      inbuf[i] = temp_sample_float_t;
    }
    filter(inbuf, out_highpass, bytes, filter_highpass);
    filter(out_highpass, out_lowpass, bytes, filter_lowpass);
    for (i = 0; i < bytes; i++) {
      out_lowpass[i] = 0.5 * out_lowpass[i];
    }
    fwrite(out_lowpass, sizeof(float), bytes, outfile);
  }
  printf("current offset %lld data_size %lld clips %i", current_offset, data_size, clips_on);
  fflush(outfile);
  free(filter_highpass);
  free(filter_lowpass);
  free(out_highpass);
  free(out_lowpass);
  free(inbuf);
  fclose(outfile);
  fclose(f);
  return 0;
}
