# GLAD Vendor Layout

Place generated GLAD OpenGL loader files here:

- `include/glad/glad.h`
- `include/KHR/khrplatform.h`
- `src/glad.c`

The `Makefile` expects `src/glad.c` and include headers under this folder.
The canonical build path uses CMake `tokamak_viewer` with `-DBUILD_VIEWER=ON`.
