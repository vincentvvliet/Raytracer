#include "texture.h"
#include <framework/image.h>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord)
{
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    const glm::ivec2 index = glm::ivec2(texCoord * glm::vec2(image.width, image.height));
    return image.pixels[index.y * image.width + (image.width - index.x) - 0.5f];
}