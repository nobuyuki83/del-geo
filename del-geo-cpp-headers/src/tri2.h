
namespace tri2 {

__device__
float area(
  const float* p0,
  const float* p1,
  const float* p2)
{
  return 0.5f * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
}


}