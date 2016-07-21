cmake_minimum_required(VERSION 3.1)

project(musal)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_CXX_STANDARD 11)

file(
	GLOB_RECURSE
	source_files
	src/*
)

include_directories(
	${CMAKE_SOURCE_DIR}/include
)

add_executable(
	musal 
	${source_files}
)

