# CUDA_ClothDynamics_PBD_GraphColoring

충돌처리 논문 : https://www.dbpia.co.kr/journal/articleDetail?nodeId=NODE11037806
 - 효율적인 BVH 빌드 방법과 R-Triangle기법으로 충돌처리 가속화

XPBD기법으로 옷감 시뮬레이션을 구현하고 GraphColoring 기법으로 가속화

# GraphColoring 빌드과정
![캡처1](https://user-images.githubusercontent.com/86860544/228162800-9cfeaa45-2ed1-44c8-8a71-c01bfbe951a8.JPG)
![캡처2](https://user-images.githubusercontent.com/86860544/228162805-86a0e1b5-5b5d-4c3b-85c8-72e909a6b05b.JPG)
![캡처3](https://user-images.githubusercontent.com/86860544/228162813-a4ad4c4c-63ca-4264-8e26-3df8e72eab8f.JPG)

# Scene
![scene](https://user-images.githubusercontent.com/86860544/228162854-b719f67a-e0f0-4d6b-ac8a-76161953c94a.gif)

# 참고문헌
 - Sean Curtis, Rasmus Tamstorf, Dinesh Manocha. Fast Collision Detection for Deformable Models using Representative-Triangles. I3D 2008, Redwood City, California, February 15–17, 2008.
 - M. Fratarcangeli1 and F. Pellacini2. Scalable Partitioning for Parallel Position Based Dynamics. EUROGRAPHICS 2015. Volume 34 (2015), Number 2
