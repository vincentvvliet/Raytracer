#include "texture.h"
#include <framework/image.h>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    if (features.enableTextureMapping) {
        int index = 3 * (texCoord[1] * image.width + texCoord[0]);
        return image.pixels[index];
    }
    
    return glm::vec3 { 0.0f, 0.0f, 0.0f };
}