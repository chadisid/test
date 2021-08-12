// compile using gcc test_pex.c -o test -lm
#include <string.h>

#include <stdio.h>

#include <stdlib.h>

#include <stdint.h>

#include <math.h>

#define PI 3.14159265359
#define SAMPLE_MIN (1 << 31)
#define SAMPLE_MAX (((unsigned)-1)>>1)
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
int check_file(const char * filename) {
  char * dot = strchr(filename, '.');
  if (!dot) {
    fprintf(stderr, "Dosent end with extension %s\n", filename);
    return -1;
  }
  if (strcmp(dot, ".raw") != 0) {
    fprintf(stderr, "Dosent end with extension raw %s\n", filename);
    return -1;
  }
  return 0;
}
int main(int argc, char ** argv) {
  int ret, num_samples = 1024, bytes, sample_rate = 48000, temp_sample, i;
  const char * outfilename, * filename;
  FILE * f, * outfile;
  int64_t data_size = 0, current_offset = 0;
  double frequency_lowpass = 300, frequency_highpass = 200, double_temp;
  size_t clips_on = 0;
  float sample_float;
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <input file> <output file>\n", argv[0]);
    exit(1);
  }
  filename = argv[1];
  ret = check_file(filename);
  if (ret < 0)
    exit(1);
  outfilename = argv[2];
  ret = check_file(outfilename);
  if (ret < 0)
    exit(1);
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
  fseek(f, 0, SEEK_END);
  data_size = ftell(f);
  fseek(f, 0, SEEK_SET);
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
  while (current_offset < data_size) {
    bytes = fread(inbuf, sizeof(float), num_samples, f);
    current_offset += sizeof(float) * bytes;
    for (i = 0; i < bytes; i++) {
      double_temp = (inbuf[i]) * (SAMPLE_MAX + 1.0);
      if (double_temp < 0) {
        if (double_temp <= SAMPLE_MIN - 0.5) {
          ++(clips_on);
          temp_sample = SAMPLE_MIN;
        } else {
          temp_sample = double_temp - 0.5;
        }
      } else {
        if (double_temp >= SAMPLE_MAX + 0.5) {
          if (double_temp > SAMPLE_MAX + 1.0) {
            ++(clips_on);
            temp_sample = SAMPLE_MAX;
          } else {
            temp_sample = SAMPLE_MAX;
          }
        } else {
          temp_sample = double_temp + 0.5;
        }
      }
      sample_float = (temp_sample) * (1.0 / (SAMPLE_MAX + 1.0));
      inbuf[i] = sample_float;
    }
    filter(inbuf, out_highpass, bytes, filter_highpass);
    filter(out_highpass, out_lowpass, bytes, filter_lowpass);
    for (i = 0; i < bytes; i++) {
      out_lowpass[i] = 0.5 * out_lowpass[i];
    }
    fwrite(out_lowpass, sizeof(float), bytes, outfile);
  }
  printf("Read %lld with file size %lld with clips %i\n", current_offset, data_size, clips_on);
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
