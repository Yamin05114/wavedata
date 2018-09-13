#include <cstdlib>
#include <iostream>
#include <stdint.h>

/*
    长度为6的列向量(其实就是一个box)拼出来的矩阵, N是box的数目
    卷积图（x，y，z）对应的数值，其实就是0->x 0->y 0->z
    那三个线段拼起来的封闭长方体里面的每一个像素之和
    Calculates sum of values given a set of 'N' coordinates stored in 'boxes'
    and a 3D integral image 'img'. 'W', 'H', and 'L' are the dimensions of the
    integral image. The coordinates should be stored in 'boxes' as N column
    vectors of size 6 (i.e. [x1, y1, z1, x2, y2, z2]). The results are stored
    in a Nx1 column vector called 'vals'.

    Arguments:
    boxes -- pointer to coordinates of boxes stored in a contiguous 6xN array
    N -- number of sets of coordinates in boxes
    img -- pointer to integral image stored in a contiguous array
    W -- width dimension of the integral image (x-axis)
    H -- height dimension of the integral image (y-axis)
    L -- length dimension of the integral image (z-axis)
    vals -- pointer to Nx1 column vector where the final values are stored
*/
extern "C" void integralImage3DVal(const uint32_t* boxes, const uint32_t N,
                                    const float* img, const uint32_t W,
                                    const uint32_t H, const uint32_t L,
                                    float* vals)
{
    uint32_t x1, y1, z1, x2, y2, z2, dx, dy, dz;
    const float* const max_ptr = img + L*W*H - 1;  // img的末尾

    for (uint32_t i = 0; i < N; ++i)  //对于每一个box
    {
        const uint32_t* box = boxes + i*6; // 指针位置正好跳六下，是下一个box
        x1 = *box - 1;            // x1 y1 z1各减一
        y1 = *(box + 1) - 1;      // x2 y2 z2保持不变
        z1 = *(box + 2) - 1;
        x2 = *(box + 3);
        y2 = *(box + 4);
        z2 = *(box + 5);

        const float* im1 = img + x1 + y1*W + z1*W*H;  
        // 这里的意思就是先增长x，在增长y，最后增长z。
        // 这个做法跟我们平时的做法是不一样的，我们在tensorflow中肯定是先增长c channel，
        // 在增长y column位置，最后增长x row位置!!!
        
        const float* im2 = img + x2 + y2*W + z2*W*H;

        // error checking to prevent accessing invalid memory
        if (im1 > max_ptr || im2 > max_ptr)
        {
            return;
        }
        else if (im1 < img || im2 < img)
        {
            return;
        }

        // get box sizes
        dx = x2 - x1;
        dy = (y2 - y1)*W;
        dz = (z2 - z1)*W*H;

        // summation based on summed area table algorithm
        vals[i] = *im2 + *(im1 + dz)
                + *(im1 + dy) + *(im1 + dx)
                - *(im2 - dy) - *(im2 - dx)
                - *(im2 - dz) - *im1;
    }
}


extern "C" void integralImage3D(float* output, const uint32_t* boxes,
                                uint32_t N, const float* image, uint32_t W,
                                uint32_t H, uint32_t L)
{
    if (boxes == NULL || image == NULL || output == NULL)
    {
        std::cout << "Invalid pointer passed into integralImage3D function." << std::endl;
        return;
    }

    integralImage3DVal(boxes, N, image, W, H, L, output);
}
