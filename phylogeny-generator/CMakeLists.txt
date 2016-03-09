project( phylogeny-generator )

cmake_minimum_required(VERSION 2.6)

ADD_DEFINITIONS(
  -std=c++11
  -w
  -lpthread
  # -ljpeg
  )

find_package(ImageMagick COMPONENTS Magick++ MagickCore)
find_package(OpenCV REQUIRED)
find_package(JPEG)

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

file(

  GLOB_RECURSE

  source_files

  src/*

  )

add_executable(

  phylogeny-generator

  ${source_files}

  )

include_directories(${ImageMagick_INCLUDE_DIRS})
include_directories(${JPEG_INCLUDE_DIRS})
target_link_libraries(phylogeny-generator ${ImageMagick_LIBRARIES} ${JPEG_LIBRARIES} ${OpenCV_LIBS} -lpthread -lX11)
