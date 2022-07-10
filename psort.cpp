#include "psort.h"
#include <omp.h>

void SequentialSort(uint32_t *data, uint32_t n);


void ParallelSort(uint32_t *data, uint32_t n, int p)
{   
    // Entry point to your sorting implementation.
    // Sorted array should be present at location pointed to by data.
    if(n<p*p) {
        SequentialSort(data, n);
        return ;
    }

    uint32_t *R = new uint32_t[p*p];
    uint32_t per_bucket = n/p;
    int extra = n%p;
    for(int i=0;i<p;i++) {
        uint32_t start = per_bucket*i;
        start = start + (i<extra ? i : extra);
        for(int j=0;j<p;j++) {
            R[i*p+j] = data[start+j];
        }
    }
    SequentialSort(R, p*p);
    uint32_t *S = new uint32_t[p-1];
    for(int i=0;i<p-1;i++){
        S[i] = R[(i+1)*p];
    }
    delete[] R;

    uint32_t *count = new uint32_t[p*p]();
    int *buck = new int[n];
    for(int i=0;i<p;i++){
        #pragma omp task firstprivate(i) 
        {
            uint32_t start = per_bucket*i;
            start = start + (i<extra ? i : extra);
            uint32_t end = start + per_bucket;
            end = end + (i<extra ? 1 : 0);
            for(uint32_t j=start;j<end;j++) {
                uint32_t val = data[j];
                for(int k=0;k<p;k++) {
                    long long l,u;
                    if(k==0) {
                        l = -1;
                        u = S[k];
                    } else if(k==p-1) {
                        l = S[k-1];
                        u = UINT32_MAX;
                    } else {
                        l = S[k-1];
                        u = S[k];
                    }
                    if(val>l && val<=u) {
                        count[k*p+i]++;
                        buck[j] = k;
                        break;
                    }
                }
            }
        }
    }
    #pragma omp taskwait
    delete[] S;

    uint32_t *pref = new uint32_t[p*(p+1)];
    for(int i=0;i<p;i++) {
        pref[i*p]=0;
        for(int j=0;j<p;j++) {
            pref[i*p+j+1]=pref[i*p+j]+count[i*p+j];
        }
    }
    
    for(int i=0;i<p;i++) {
        for(int j=1;j<p;j++) {
            count[i*p]+=count[i*p+j];
            count[i*p+j]=0;
        }
    }

    uint32_t *precn = new uint32_t[p+1];
    precn[0]=0;
    for(int i=0;i<p;i++) {
        precn[i+1]=precn[i]+count[i*p];
        count[i*p]=0;
    }

    uint32_t *dum = new uint32_t[n];
    for(int i=0;i<p;i++) {
        #pragma omp task firstprivate(i) 
        {
            uint32_t start = per_bucket*i;
            start = start + (i<extra ? i : extra);
            uint32_t end = start + per_bucket;
            end = end + (i<extra ? 1 : 0);
            for(uint32_t j=start;j<end;j++) {
                dum[precn[buck[j]]+pref[buck[j]*p+i]+count[buck[j]*p+i]] = data[j];
                count[buck[j]*p+i]++;
            }
        }
    }
    #pragma omp taskwait
    delete[] buck;
    delete[] pref;
    delete[] count;
    for(int i=0;i<p;i++) {
        #pragma omp task firstprivate(i) 
        {
            uint32_t start = per_bucket*i;
            start = start + (i<extra ? i : extra);
            uint32_t end = start + per_bucket;
            end = end + (i<extra ? 1 : 0);
            for(uint32_t j=start;j<end;j++) {
                data[j]=dum[j];
            }
        }
    }
    #pragma omp taskwait
    delete[] dum;

    uint32_t threshold = (2*(uint32_t)(n))/p;
    for(int i=0;i<p;i++){
        #pragma omp task firstprivate(i)
        {
            uint32_t cn = precn[i+1]-precn[i];
            if(cn<threshold) {
                SequentialSort(data+precn[i],cn);
            } else {
                ParallelSort(data+precn[i],cn,p);
            }
        }
    }
    #pragma omp taskwait
    delete[] precn;

    // uint32_t itr=0;
    // while(itr<n) {
    //     if(buck[itr]!=p) {
    //         uint32_t dst=precn[buck[itr]]+count[buck[itr]*p];
    //         if(itr==dst) {
    //             count[buck[itr]*p]++;
    //             buck[itr]=p;
    //             itr++;
    //         } else {
    //             uint32_t tmp=data[dst];
    //             data[dst]=data[itr];
    //             data[itr]=tmp;
    //             count[buck[itr]*p]++;
    //             buck[itr]=buck[dst];
    //             buck[dst]=p;
    //         }
    //     } else {
    //         itr++;
    //     }
    // }
    // delete[] buck;

    // uint32_t threshold = (2*(uint32_t)(n))/p;
    // for(int i=0;i<p;i++){
    //     #pragma omp task firstprivate(i)
    //     {
    //         if(count[i*p]<threshold) {
    //             SequentialSort(data+precn[i],count[i*p]);
    //         } else {
    //             ParallelSort(data+precn[i],count[i*p],p);
    //         }
    //     }
    // }
    // #pragma omp taskwait
    // delete[] precn;
    // delete[] count;

}

void SequentialSort(uint32_t *data, uint32_t n) {
    if(n<=1) {
        return;
    }
    uint32_t mid = n/2;

    SequentialSort(data,mid);
    SequentialSort(data+mid,n-mid);

    uint32_t *tmp = new uint32_t[n];
    uint32_t i=0,j=mid,k=0;
    while(i<mid && j<n) {
        if(data[i]<=data[j]) {
            tmp[k++] = data[i++];
        } else {
            tmp[k++] = data[j++];
        }
    }
    while(i<mid) {
        tmp[k++] = data[i++];
    }
    while(j<n) {
        tmp[k++] = data[j++];
    }
    k=0;
    while(k<n) {
        data[k] = tmp[k];
        k++;
    }
    delete[] tmp;
}




// void ParallelSort2(uint32_t *data, uint32_t n, int p)
// {   
//     // Entry point to your sorting implementation.
//     // Sorted array should be present at location pointed to by data.
//     if(n<p*p) {
//         SequentialSort(data, n);
//         return ;
//     }

//     uint32_t *R = new uint32_t[p*p];
//     uint32_t per_bucket = n/p;
//     int extra = n%p;
//     for(int i=0;i<p;i++) {
//         uint32_t start = per_bucket*i;
//         start = start + (i<extra ? i : extra);
//         for(int j=0;j<p;j++) {
//             R[i*p+j] = data[start+j];
//         }
//     }
//     SequentialSort(R, p*p);
//     uint32_t *S = new uint32_t[p-1];
//     for(int i=0;i<p-1;i++){
//         S[i] = R[(i+1)*p];
//     }
//     delete[] R;
//     uint32_t threshold = (2*(uint32_t)(n))/p;
//     uint32_t *count = new uint32_t[p]();
//     int *buck = new int[n]();
//     uint32_t **newarr = new uint32_t*[n];
//     for(int i=0;i<p;i++){
//         #pragma omp task firstprivate(i) 
//         {
//             long long l,u;
//             if(i==0) {
//                 l = -1;
//                 u = S[i];
//             } else if(i==p-1) {
//                 l = S[i-1];
//                 u = UINT32_MAX;
//             } else {
//                 l = S[i-1];
//                 u = S[i];
//             }
//             for(uint32_t j=0;j<n;j++) {
//                 uint32_t val = data[j];
//                 if(val>l && val<=u) {
//                     count[i]++;
//                     buck[j] = i+1;
//                 }
//             }
//             newarr[i] = new uint32_t[count[i]];
//             uint32_t k=0;
//             for(uint32_t j=0;j<n;j++) {
//                 if(buck[j]-1==i) {
//                     newarr[i][k] = data[j];
//                     k++;
//                 }
//             }
//             if(count[i]<threshold) {
//                 SequentialSort(newarr[i],count[i]);
//             } else {
//                 ParallelSort2(newarr[i],count[i],p);
//             }
//         }
//     }
//     #pragma omp taskwait
//     delete[] S;
//     delete[] buck;

//     uint32_t total=0;
//     for(int i=0;i<p;i++){
//         #pragma omp task firstprivate(i, total)
//         {
//             for(uint32_t j=0;j<count[i];j++){
//                 data[total+j] = newarr[i][j];
//             }
//             delete[] newarr[i];
//         }
//         total+=count[i];
//     }
//     #pragma omp taskwait
//     delete[] newarr;
//     delete[] count;

// }