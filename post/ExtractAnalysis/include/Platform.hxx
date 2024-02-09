/**
 * File:   Platform.hxx
 * Author: akirby
 *
 * Created on September 29, 2021, 2:50 PM
 */

#ifndef PLATFORM_HXX
#define PLATFORM_HXX

#include <mpi.h>
#include <occa.hpp>
#include <unistd.h>
#include <cstring>

#ifdef __DG_CUDA__
#  include <cuda_runtime.h>
#endif

#ifdef __DG_HIP__
#  include <hip/hip_runtime.h>
#endif

#define HIP_MODE 0
#define CUDA_MODE 1
#define OpenCL_MODE 2
#define OpenMP_MODE 3
#define Serial_MODE 4

class Platform {
  public:
    MPI_Comm comm;

    occa::properties props;
    occa::device device;

    int rank,nrank;
    int thread_model;

    /* constructors */
    Platform(MPI_Comm _comm,int _thread_model,int device_id=0,int plat=0):
        comm(_comm),thread_model(_thread_model)
    {
        MPI_Comm_rank(_comm, &rank);
        MPI_Comm_size(_comm, &nrank);

        if(rank==0){
            fprintf(stdout,"Thread Model: %s, device_id: %d",
                    thread_model == HIP_MODE    ? "HIP":
                    thread_model == CUDA_MODE   ? "CUDA":
                    thread_model == OpenCL_MODE ? "OpenCL":
                    thread_model == OpenMP_MODE ? "OpenMP":"Serial",
                    device_id);
            (thread_model == OpenCL_MODE) ? printf(", platform: %d\n",plat):printf("\n");
        }
        DeviceConfig(thread_model,device_id,plat);
    }
   ~Platform(){}

    /* methods */
    occa::kernel buildKernel(const std::string &fileName,
                             const std::string &kernelName,
                             const occa::json  &props = occa::json()){
        occa::kernel kernel;

        /* build on root first */
        if(!rank) kernel = device.buildKernel(fileName,kernelName,props);
        MPI_Barrier(comm);

        /* remaining ranks find the cached version (ideally) */
        if(rank) kernel = device.buildKernel(fileName,kernelName,props);
        MPI_Barrier(comm);

        return kernel;
    }

    template <class T>
    occa::memory malloc(const occa::dim_t entries,
                        const void *src = NULL,
                        const occa::json &prop = occa::properties()) {
        return device.malloc<T>(entries, src, prop);
    }

    template <class T>
    occa::memory malloc(const occa::dim_t entries,
                        const occa::memory &src,
                        const occa::json &prop = occa::properties()) {
        return device.malloc<T>(entries, src, prop);
    }

    template <class T>
    occa::memory malloc(const occa::dim_t entries,
                        const occa::json &prop) {
        return device.malloc<T>(entries, prop);
    }

    void *hostMalloc(const size_t bytes,
                     const void *src,
                     occa::memory &h_mem){
        occa::properties hostProp;
        hostProp["host"] = true;
        h_mem = device.malloc(bytes, src, hostProp);
        return h_mem.ptr();
    }

    void memcpyToSymbol(const void *dst,
                        const void *src,
                        const size_t bytes){
#ifdef __DG_HIP__
        if(thread_model == HIP_MODE) hipMemcpyToSymbol(dst, src, bytes);
#endif
#ifdef __DG_CUDA__
        if(thread_model == CUDA_MODE) cudaMemcpyToSymbol(dst, src, bytes);
#endif
    }

  private:
    void DeviceConfig(int thread_model,int device_id,int plat);
};

#endif /* PLATFORM_HXX */