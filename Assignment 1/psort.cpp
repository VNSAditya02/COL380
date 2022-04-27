#include "psort.h"
#include <omp.h>

uint32_t* SequentialSort(uint32_t *v, uint32_t l, uint32_t r){
    uint32_t *result;
    if(l > r){
        return result;
    }
    result = new uint32_t[r - l + 1];
    if(l == r){
        result[0] = v[l];
        return result;
    }
    else if(r - l == 1){
        if(v[l] > v[r]){
            result[0] = v[r];
            result[1] = v[l];
        }
        else{
            result[0] = v[l];
            result[1] = v[r];
        }
        return result;
    }

    uint32_t *result1 = SequentialSort(v, l, l + (r - l)/2);
    uint32_t *result2 = SequentialSort(v, 1 + l + (r - l)/2, r);
    uint32_t n1 = 1 + (r - l)/2, n2 = r - (l + (r - l)/2);
    uint32_t lp = 0, rp = 0;
    uint32_t k = 0;
    while (lp < n1 && rp < n2){
        if (result1[lp] <= result2[rp]){
            result[k] = result1[lp];
            lp++;
            k++;
        }
        else{
            result[k] = result2[rp];
            rp++;
            k++;
        }
    }
    while (lp < n1) {
        result[k] = result1[lp];
        lp++;
        k++;
    }

    while (rp < n2) {
        result[k] = result2[rp];
        rp++;
        k++;
    }
    return result;
}

void countV(uint32_t *data, uint32_t *splitters, uint32_t *count[], uint32_t start, uint32_t end, int p, int id)
{
    for(uint32_t i = start; i < end; i++){
        uint32_t key = data[i];
        int pos = -1;
        int low = 0, high = p - 2;
        while (low <= high) {
            int mid = low + (high - low + 1) / 2;
            uint32_t midVal = splitters[mid];
            if (midVal < key) {
                low = mid + 1;
            }
            else if (midVal > key) {
                pos = mid;
                high = mid - 1;
            }
            else if (midVal == key) {
                low = mid + 1;
            }
        }
        if(pos == -1){
            pos = p - 1;
        }
        count[pos][id]++;
    }
}

void workOntask(uint32_t *data, uint32_t *B[][25], uint32_t *splitters, uint32_t *count[], uint32_t start, uint32_t end, int p, int id)
{
    for(uint32_t i = start; i < end; i++){
        uint32_t key = data[i];
        int pos = -1;
        int low = 0, high = p - 2;
        while (low <= high) {
            int mid = low + (high - low + 1) / 2;
            uint32_t midVal = splitters[mid];
            if (midVal < key) {
                low = mid + 1;
            }
            else if (midVal > key) {
                pos = mid;
                high = mid - 1;
            }
            else if (midVal == key) {
                low = mid + 1;
            }
        }
        if(pos == -1){
            pos = p - 1;
        }
        B[pos][id][count[pos][id]] = data[i];
        count[pos][id]++;
    }
}

void joinAll(uint32_t *data, uint32_t *final_B[], uint32_t *final_count, uint32_t *B[][25], uint32_t *count[], int i, int parts)
{
    for(int j = 0; j < parts; j++){
        for(uint32_t k = 0; k < count[i][j]; k++){
            final_B[i][final_count[i]] = B[i][j][k];
            final_count[i]++;
        }
    }
}

void concatenate(uint32_t *data, uint32_t *final_B[], uint32_t ni, uint32_t start, int i)
{
    uint32_t k = start;
    for(uint32_t j = 0; j < ni; j++){
        data[k] = final_B[i][j];
        k++;
    }
}

void sortAgain(uint32_t *final_B[], uint32_t ni, uint32_t n, int p, int i)
{
    if(ni < 2*n/p){
        uint32_t *result = SequentialSort(final_B[i], 0, ni - 1);
        for(uint32_t j = 0; j < ni; j++){
            final_B[i][j] = result[j];
        }
    }
    else{
        ParallelSort(final_B[i], ni, p);
    }
}

void ParallelSort(uint32_t *data, uint32_t n, int p)
{
    // Entry point to your sorting implementation.
    // Sorted array should be present at location pointed to by data.

    // Base Case
    if(n <= p*p){
        uint32_t *R = SequentialSort(data, 0, n - 1);
        for(uint32_t i = 0; i < n; i++){
            data[i] = R[i];
        }
        return;
    }

    // Step 1
    uint32_t *v = new uint32_t[p*p];
    uint32_t s = 0, k = 0;
    for(int i = 0; i < p; i++){
    	for(uint32_t j = s; j < s + p; j++){
            v[k] = data[j];
            k++;
        }
        if(i < n % p){
            s += (n / p) + 1;
        }
        else{
            s += (n / p);
        }
        
    }

    // Step 2
    uint32_t *R = SequentialSort(v, 0, p*p - 1);
    
    // Step 3
    uint32_t *splitters = new uint32_t[p - 1];
    for(int i = 0; i < p - 1; i++){
        splitters[i] = R[(i + 1)*p];
    }
    
    // Step 4
    int parts = 25; // Change parts in different places
    uint32_t *tmp_count[p];
    for(int i = 0; i < p; i++){
        tmp_count[i] = new uint32_t[25];
        for(int j = 0; j < parts; j++){
            tmp_count[i][j] = 0;
        }
    }

    s = 0;
    uint32_t e;
    for(int i = 0; i < parts; i++){
        if(i < n % parts){
            e = s + (n / parts) + 1;
        }
        else{
            e = s + (n / parts);
        }
        #pragma omp task
        {
            countV(data, splitters, tmp_count, s, e, p, i);
        }
        s = e;
    }
    #pragma omp taskwait

    uint32_t *count[p];
    uint32_t *total_tmp_count = new uint32_t[p];
    uint32_t *B[p][25];

    for(int i = 0; i < p; i++){
        total_tmp_count[i] = 0;
    }

    for(int i = 0; i < p; i++){
        count[i] = new uint32_t[25];
        for(int j = 0; j < parts; j++){
            B[i][j] = new uint32_t[tmp_count[i][j]];
            count[i][j] = 0;
            total_tmp_count[i] += tmp_count[i][j];
        }
    }
    
    s = 0;
    for(int i = 0; i < parts; i++){
        if(i < n % parts){
            e = s + (n / parts) + 1;
        }
        else{
            e = s + (n / parts);
        }
        #pragma omp task
        {
            workOntask(data, B, splitters, count, s, e, p, i);
        }
        s = e;
    }
    #pragma omp taskwait

    // Step 5
    uint32_t *final_count = new uint32_t[p];
    uint32_t *final_B[p];
    for(int i = 0; i < p; i++){
        final_count[i] = 0;
        final_B[i] = new uint32_t[total_tmp_count[i]];
        #pragma omp task
        {
            joinAll(data, final_B, final_count, B, count, i, parts);
        }
    }
    #pragma omp taskwait

    for(int i = 0; i < p; i++){
        uint32_t ni = final_count[i];
        #pragma omp task
        {
            sortAgain(final_B, ni, n, p, i);
        }
    }
    #pragma omp taskwait

    uint32_t *range_count = new uint32_t[p];
    range_count[0] = 0;
    for(int i = 1; i < p; i++){
        range_count[i] = range_count[i - 1] + final_count[i - 1];
    }

   // Step 6
    k = 0;
    for(int i = 0; i < p; i++){
        uint32_t ni = final_count[i];
        uint32_t start = range_count[i];
        #pragma omp task
        {
            concatenate(data, final_B, ni, start, i);
        }
    }
}
