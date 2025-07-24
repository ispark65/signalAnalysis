void detect_automatic_reflection_points(double* signal, int size, double sample_interval_s, int* out_indices, int* out_count, int max_points);

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_SAMPLES 10000
#define SMOOTH_WINDOW 10
#define MAX_LINE_LENGTH 1024

double rf_data[MAX_SAMPLES];
double envelope[MAX_SAMPLES];

double smooth(double* arr, int idx, int size, int window) {
    double sum = 0.0;
    int count = 0;
    for (int i = idx - window/2; i <= idx + window/2; i++) {
        if (i >= 0 && i < size) {
            sum += arr[i];
            count++;
        }
    }
    return (count > 0) ? (sum / count) : 0.0;
}

void generate_output_filename(const char* input, char* output, size_t maxlen, const char* new_ext) {
    strncpy(output, input, maxlen - 1);
    output[maxlen - 1] = '\0';

    char *dot = strrchr(output, '.');
    if (dot) {
        *dot = '\0';
    }
    strncat(output, "-chg", maxlen - strlen(output) - 1);
    strncat(output, new_ext, maxlen - strlen(output) - 1);
}

int is_csv_file(const char* filename) {
    const char *ext = strrchr(filename, '.');
    return ext && strcmp(ext, ".csv") == 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <input_file> [output_file]\n", argv[0]);
        return 1;
    }

    const char *input_filename = argv[1];
    char output_filename[256];

    if (argc >= 3) {
        strncpy(output_filename, argv[2], sizeof(output_filename) - 1);
        output_filename[sizeof(output_filename) - 1] = '\0';
    } else {
        const char* ext = is_csv_file(input_filename) ? ".csv" : ".txt";
        generate_output_filename(input_filename, output_filename, sizeof(output_filename), ext);
    }

    int count = 0;
    if (is_csv_file(input_filename)) {
        FILE *fp = fopen(input_filename, "r");
        if (!fp) {
            perror("Failed to open CSV input file");
            return 1;
        }

        char line[MAX_LINE_LENGTH];
        int line_number = 0;

        while (fgets(line, sizeof(line), fp)) {
            line_number++;
            if (line_number < 6) continue;  // Skip first 5 lines

            char *token = strtok(line, ",");
            if (!token) continue;
            token = strtok(NULL, ",");  // Second column
            if (!token) continue;

            rf_data[count] = atof(token);
            count++;
            if (count >= MAX_SAMPLES) break;
        }

        fclose(fp);
    } else {
        FILE *fp = fopen(input_filename, "r");
        if (!fp) {
            perror("Failed to open TXT input file");
            return 1;
        }

        while (count < MAX_SAMPLES && fscanf(fp, "%lf", &rf_data[count]) == 1) {
            count++;
        }

        fclose(fp);
    }

    for (int i = 0; i < count; i++) {
        envelope[i] = fabs(rf_data[i]);
    }

    double smoothed[MAX_SAMPLES];
    for (int i = 0; i < count; i++) {
        smoothed[i] = smooth(envelope, i, count, SMOOTH_WINDOW);
    }

    FILE *out = fopen(output_filename, "w");
    if (!out) {
        perror("Failed to create output file");
        return 1;
    }

    if (is_csv_file(output_filename)) {
        fprintf(out, "Index,SmoothedEnvelope\n");
        for (int i = 0; i < count; i++) {
            fprintf(out, "%d,%f\n", i, smoothed[i]);
        }
    } else {
        for (int i = 0; i < count; i++) {
            fprintf(out, "%d\t%f\n", i, smoothed[i]);
        }
    }

    fclose(out);
    printf("Conversion completed: %s\n", output_filename);

    int reflection_indices[10];
    int reflection_count = 0;
    double sample_interval_s = 10e-9;

    detect_automatic_reflection_points(smoothed, count, sample_interval_s, reflection_indices, &reflection_count, 10);

    printf("Automatically detected reflection points based on gradient rise and peak separation:\n");
    for (int i = 0; i < reflection_count; i++) {
        printf("Reflection %d: Index = %d\n", i + 1, reflection_indices[i]);
    }

    return 0;

}


#include <float.h>
#include <stdlib.h>
#include <math.h>

// 자동 반사 신호 탐지 (첫 반사 이후 3μs, 이후 중복 방지 위해 최소 30샘플 간격 적용)
void detect_automatic_reflection_points(double* signal, int size, double sample_interval_s, int* out_indices, int* out_count, int max_points) {
    int last_peak_index = -100000;
    int detected = 0;
    double threshold = 0.1;
    int min_interval_samples_first = (int)(3e-6 / sample_interval_s);  // 3μs
    int min_interval_samples_repeat = 30;  // 이후 중복 방지를 위한 최소 이격 샘플 수

    for (int i = 1; i < size - 1 && detected < max_points; i++) {
        if (signal[i] > signal[i-1] && signal[i] > signal[i+1] && signal[i] > threshold) {
            int min_interval = (detected == 0) ? min_interval_samples_first : min_interval_samples_repeat;
            if (i - last_peak_index >= min_interval) {
                // 앞쪽으로 기울기 최대 지점 찾기
                int start = (i - 30 > 0) ? i - 30 : 0;
                double max_grad = -DBL_MAX;
                int max_grad_idx = start;
                for (int j = start + 1; j < i; j++) {
                    double grad = signal[j] - signal[j - 1];
                    if (grad > max_grad) {
                        max_grad = grad;
                        max_grad_idx = j;
                    }
                }
                // 중복 방지: 이전에 이미 선택된 인덱스와 유사하면 제외
                int is_duplicate = 0;
                for (int k = 0; k < detected; k++) {
                    if (abs(out_indices[k] - max_grad_idx) < min_interval_samples_repeat) {
                        is_duplicate = 1;
                        break;
                    }
                }
                if (!is_duplicate) {
                    out_indices[detected++] = max_grad_idx;
                    last_peak_index = i;
                }
            }
        }
    }
    *out_count = detected;
}
