#include "helpers.h"
#include <math.h>

// Convert image to grayscale
void grayscale(int height, int width, RGBTRIPLE image[height][width])
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            // Calculate the average of the RGB values of each pixel
            RGBTRIPLE pixel = image[i][j];
            BYTE average_rgb = (BYTE)round((pixel.rgbtBlue + pixel.rgbtGreen + pixel.rgbtRed) / 3.0);
            // Set each RGB value to the calculated average
            image[i][j].rgbtBlue = average_rgb;
            image[i][j].rgbtGreen = average_rgb;
            image[i][j].rgbtRed = average_rgb;
        }
    }
    return;
}

// Reflect image horizontally
void reflect(int height, int width, RGBTRIPLE image[height][width])
{
    for (int i = 0; i < height; i++)
    {
        // Swap each pixel in image with its corresponding pixel from the opposite edge, inverting the image horizontally
        for (int j = 0; j < width / 2; j++)
        {
            RGBTRIPLE temp_pixel = image[i][j];
            image[i][j] = image[i][width - 1 - j];
            image[i][width - 1 - j] = temp_pixel;
        }
    }
    return;
}

// Blur image
void blur(int height, int width, RGBTRIPLE image[height][width])
{
    // Create a copy of image array called temp to ease actual process below
    RGBTRIPLE temp[height][width];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            temp[i][j] = image[i][j];
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            // For each pixel in temp, check its 3x3 surroundings and store the sum RGB values of existing neighbor pixels
            int red_sum = 0, green_sum = 0, blue_sum = 0;
            float count = 0.0;
            for (int k = i - 1; k < i + 2; k++)
            {
                if (k >= 0 && k < height)
                {
                    for (int l = j - 1; l < j + 2; l++)
                    {
                        if (l >= 0 && l < width)
                        {
                            red_sum += temp[k][l].rgbtRed;
                            green_sum += temp[k][l].rgbtGreen;
                            blue_sum += temp[k][l].rgbtBlue;
                            count++;
                        }
                    }
                }
            }
            // Calculate the average of the total surrounding values for each RGB value and assign it to the central pixel in image
            image[i][j].rgbtRed = (BYTE)round(red_sum / count);
            image[i][j].rgbtGreen = (BYTE)round(green_sum / count);
            image[i][j].rgbtBlue = (BYTE)round(blue_sum / count);
        }
    }

    return;
}

// Detect edges
void edges(int height, int width, RGBTRIPLE image[height][width])
{
    // Create a copy of image array called temp to ease actual process below
    RGBTRIPLE temp[height][width];
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            temp[i][j] = image[i][j];
        }
    }

    // Create the Gx kernel array
    int Gx[3][3];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Gx[i][j] = ((i % 2) + 1) * (j - 1);
        }
    }

    // Create the Gy kernel array
    int Gy[3][3];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Gy[i][j] = ((j % 2) + 1) * (i - 1);
        }
    }

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            // For each pixel in temp, check its 3x3 surroundings and store the kernel values of existing neighbor pixels
            int red_gx = 0, green_gx = 0, blue_gx = 0;
            int red_gy = 0, green_gy = 0, blue_gy = 0;
            for (int k = i - 1; k < i + 2; k++)
            {
                if (k >= 0 && k < height)
                {
                    for (int l = j - 1; l < j + 2; l++)
                    {
                        if (l >= 0 && l < width)
                        {
                            // We multiply each RGB value value with its positionally corresponding Gx and Gy entries
                            // and keep track of the sum
                            red_gx += temp[k][l].rgbtRed * Gx[k - (i - 1)][l - (j - 1)];
                            red_gy += temp[k][l].rgbtRed * Gy[k - (i - 1)][l - (j - 1)];
                            green_gx += temp[k][l].rgbtGreen * Gx[k - (i - 1)][l - (j - 1)];
                            green_gy += temp[k][l].rgbtGreen * Gy[k - (i - 1)][l - (j - 1)];
                            blue_gx += temp[k][l].rgbtBlue * Gx[k - (i - 1)][l - (j - 1)];
                            blue_gy += temp[k][l].rgbtBlue * Gy[k - (i - 1)][l - (j - 1)];
                        }
                    }
                }
            }
            // Calculate the Sobel Operator for each RGB value
            float redEdge = sqrt(red_gx * red_gx + red_gy * red_gy);
            float greenEdge = sqrt(green_gx * green_gx + green_gy * green_gy);
            float blueEdge = sqrt(blue_gx * blue_gx + blue_gy * blue_gy);

            // Assign the Sobel Operator (maxed out at 255) to the RGB value of the current pixel in image
            image[i][j].rgbtRed = (BYTE)round(fmin(redEdge, 255.0));
            image[i][j].rgbtGreen = (BYTE)round(fmin(greenEdge, 255.0));
            image[i][j].rgbtBlue = (BYTE)round(fmin(blueEdge, 255.0));
        }
    }
    return;
}
