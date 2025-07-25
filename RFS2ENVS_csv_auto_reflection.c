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

// 노이즈 추정을 위한 평균 및 표준편차 계산 헬퍼 함수
void calculate_stats(double* signal, int size, double* mean, double* std_dev) {
    if (size <= 1) { // 표준편차를 계산하려면 최소 2개의 포인트가 필요
        *mean = (size == 1) ? signal[0] : 0;
        *std_dev = 0;
        return;
    }

    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += signal[i];
    }
    *mean = sum / size;

    double sum_sq_diff = 0.0;
    for (int i = 0; i < size; i++) {
        sum_sq_diff += (signal[i] - *mean) * (signal[i] - *mean);
    }
    // 표본 표준편차(N-1로 나눔)를 사용하여 더 나은 추정치 계산
    *std_dev = sqrt(sum_sq_diff / (size - 1));
}


// 최종 하이브리드 동적 분석 알고리즘
void detect_automatic_reflection_points(double* signal, int size, double sample_interval_s, int* out_indices, int* out_count, int max_points) {
    *out_count = 0;
    if (size < 50) return;

    // --- 1. 초기 설정 및 노이즈 추정 ---
    int temp_indices[max_points];
    int temp_count = 0;

    int noise_sample_size = size / 10 > 100 ? 100 : (size / 10 < 20 ? 20 : size / 10);
    if (noise_sample_size >= size) noise_sample_size = size - 1;

    double noise_mean, noise_std_dev;
    calculate_stats(signal, noise_sample_size, &noise_mean, &noise_std_dev);

    // --- 2. 파라미터 정의 ---
    const double K_SENSITIVITY = 5.0;
    const double SUSTAIN_FACTOR = 0.5;      // -6dB, 초기 구간 블라인드 존 설정용
    const int MIN_PEAK_INTERVAL = 30;
    const double LATE_ZONE_THRESHOLD_RATIO = 0.3; // 최대 진폭의 30% 이하부터 후기 구간으로 간주
    const double REJECTION_RATIO = 0.15;    // 최종 필터링 비율 (최대 진폭의 15%)

    double adaptive_threshold = noise_mean + K_SENSITIVITY * noise_std_dev;
    if (adaptive_threshold < 0.05) adaptive_threshold = 0.05;

    // --- 3. 하이브리드 탐색 루프 ---
    int current_search_idx = 1;
    double max_amplitude_so_far = 0.0;
    int analysis_mode = 0; // 0: 초기 구간(고해상도), 1: 후기 구간(고안정성)

    while (current_search_idx < size - 1 && temp_count < max_points) {
        int cross_idx = -1;
        for (int i = current_search_idx; i < size; i++) {
            if (signal[i] > adaptive_threshold) {
                cross_idx = i;
                break;
            }
        }
        if (cross_idx == -1) break;

        int search_window = 50;
        int search_end = cross_idx + search_window;
        if (search_end >= size) search_end = size - 1;

        double max_grad = -DBL_MAX;
        int leading_edge_idx = cross_idx;
        for (int i = cross_idx; i <= search_end; i++) {
            if (i > 0) {
                double grad = signal[i] - signal[i - 1];
                if (grad > max_grad) { max_grad = grad; leading_edge_idx = i; }
            }
        }

        int peak_idx = leading_edge_idx;
        for (int i = leading_edge_idx; i <= search_end; i++) {
            if (signal[i] > signal[peak_idx]) { peak_idx = i; }
        }

        if (temp_count > 0 && (leading_edge_idx - temp_indices[temp_count - 1]) < MIN_PEAK_INTERVAL) {
            current_search_idx = peak_idx + 1;
            continue;
        }

        double current_amplitude = signal[peak_idx];
        if (current_amplitude > max_amplitude_so_far) {
            max_amplitude_so_far = current_amplitude;
        }

        // 모드 전환 체크: 현재 피크가 최대 피크의 30%보다 작으면 후기 구간 모드로 전환
        if (analysis_mode == 0 && current_amplitude < max_amplitude_so_far * LATE_ZONE_THRESHOLD_RATIO) {
            analysis_mode = 1; // 후기 구간(고안정성) 모드로 전환
        }

        temp_indices[temp_count++] = leading_edge_idx;

        // 다음 탐색 위치 결정 (하이브리드 로직)
        if (analysis_mode == 0) { // 초기 구간: -6dB 감쇠 지점까지 블라인드
            double sustain_threshold = current_amplitude * SUSTAIN_FACTOR;
            int echo_end_idx = peak_idx;
            for (int i = peak_idx + 1; i < size; i++) {
                if (signal[i] < sustain_threshold) { echo_end_idx = i; break; }
                if (i == size - 1) { echo_end_idx = i; }
            }
            current_search_idx = echo_end_idx + 1;
        } else { // 후기 구간: 고정된 최소 간격만큼만 블라인드
            current_search_idx = peak_idx + MIN_PEAK_INTERVAL;
        }
    }

    // --- 4. 최종 상대 진폭 필터링 ---
    if (temp_count > 1) {
        double rejection_threshold = max_amplitude_so_far * REJECTION_RATIO;
        
        // 첫 번째 기준점은 항상 유지
        out_indices[(*out_count)++] = temp_indices[0];

        for (int i = 1; i < temp_count; i++) {
            // 리딩 엣지에서의 신호 크기를 기준으로 필터링
            if (signal[temp_indices[i]] >= rejection_threshold) {
                out_indices[(*out_count)++] = temp_indices[i];
            }
        }
    } else if (temp_count == 1) { // 반사점이 하나만 있으면 그대로 사용
        out_indices[0] = temp_indices[0];
        *out_count = 1;
    }
}
