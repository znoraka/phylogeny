project( tree-extractor )

cmake_minimum_required(VERSION 2.6)

ADD_DEFINITIONS(
  -std=c++11
  -w
  -lpthread
  -Ofast
  # -g
  )

find_package(ImageMagick COMPONENTS Magick++ MagickCore MagickWand)
find_package(JPEG)

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

file(

  GLOB_RECURSE

  source_files

  src/*

  )

add_executable(

  tree-extractor

  ${source_files}

  )

include_directories(${ImageMagick_INCLUDE_DIRS})
include_directories(${JPEG_INCLUDE_DIRS})
target_link_libraries(tree-extractor ${ImageMagick_LIBRARIES} ${JPEG_LIBRARIES} -lpthread -lX11)
