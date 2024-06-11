cmake_minimum_required(VERSION 3.18) 

project(shamir VERSION 1.0.0)

#version standard
Set(CMAKE_C_STANDARD 11)
Set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#print some message  eg. variable 
message(CMAKE_CXX_STANDARD="${CMAKE_CXX_STANDARD}")

set(MY_DEV_LIB_PATH "D:/workspace/dev")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(ARCH_TYPE "${CMAKE_SYSTEM_PROCESSOR}")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(ARCH_TYPE "${CMAKE_SYSTEM_PROCESSOR}")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Android")
  set(ARCH_TYPE "${CMAKE_ANDROID_ARCH_ABI}")
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  message("64bit")
  set(BUILD_BIT "x64")
else()
  # 暂时用不到32位
  message("32bit")
  set(BUILD_BIT "x86")
endif()


#static lib output path  windows: .a .lib   unix: .a
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${ARCH_TYPE}/$<IF:$<CONFIG:Debug>,Debug,Release>/lib")
# 目标链接库的存放位置 unix: .so
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${ARCH_TYPE}/$<IF:$<CONFIG:Debug>,Debug,Release>/lib")
# windows: .exe .dll   unix: 可执行文件
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${ARCH_TYPE}/$<IF:$<CONFIG:Debug>,Debug,Release>/bin")

# 通用宏定义
list(APPEND EXTRA_DEFINES)
# 编译选项 目前只支持gcc和clang, msvc用CMAKE_CXX_FLAGS
list(APPEND EXTRA_COMPILE_FLAGS)
# 依赖的库的路径
list(APPEND EXTRA_LIB_PATH)
# 依赖的库
list(APPEND EXTRA_LIB)

list(APPEND EXTRA_INCLUDE  
                           "${PROJECT_SOURCE_DIR}/include" # 本工程的头文件
                           "${MY_DEV_LIB_PATH}/include" # 第三方库的头文件
                           )

list(APPEND EXTRA_SRC "${PROJECT_SOURCE_DIR}/src/shamir.cpp")

option(DLL "build dll" ON)

#link lib 后面应该完善，操作系统和编译器组合判断
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  message("windows")
# 根据编译器的不同，设置不同的编译选项和依赖库
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message("mingw")
    # 方便生成compile_commands.json
    list(APPEND EXTRA_INCLUDE "D:/workspace/dev/vcpkg/installed/${BUILD_BIT}-windows/include")
    list(APPEND EXTRA_LIB_PATH "${MY_DEV_LIB_PATH}/unixlib")
    list(APPEND EXTRA_COMPILE_FLAGS "-Wall")

    list(APPEND EXTRA_LIB gtest pthread)

  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("clang")

  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    message("msvc")
    # set(CMAKE_TOOLCHAIN_FILE "${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
    set(CMAKE_TOOLCHAIN_FILE "D:/workspace/dev/vcpkg/scripts/buildsystems/vcpkg.cmake")
    set(CMAKE_PREFIX_PATH "D:/workspace/dev/vcpkg/installed/${BUILD_BIT}-windows")
    set(VCPKG_FEATURE_FLAGS "versions")

    # 第三方库的头文件（vcpkg安装的）
    if(DLL)
        list(APPEND EXTRA_INCLUDE "${CMAKE_PREFIX_PATH}/include")
    else()
        list(APPEND EXTRA_INCLUDE "${CMAKE_PREFIX_PATH}-static/include")
    endif()

    # 第三方库的路径（vcpkg安装的）
    # list(APPEND EXTRA_LIB_PATH "${CMAKE_PREFIX_PATH}/$<$<CONFIG:Debug>:debug/>bin")
    # 区分debug和release, 多配置编译器不能直接使用CMAKE_BUILD_TYPE
    list(APPEND EXTRA_LIB_PATH "${CMAKE_PREFIX_PATH}-static/$<$<CONFIG:Debug>:debug/>lib")
    list(APPEND EXTRA_LIB_PATH "${CMAKE_PREFIX_PATH}/$<$<CONFIG:Debug>:debug/>lib")
    # 自己编译的库的路径(msvc编译器)
    list(APPEND EXTRA_LIB_PATH "${MY_DEV_LIB_PATH}/msvclib")

    # list(APPEND EXTRA_DEFINES "_LOG_")
    # list(APPEND EXTRA_COMPILE_FLAGS "$<IF:$<CONFIG:Debug>,MTd,MT>")
    list(APPEND EXTRA_COMPILE_FLAGS "/MD$<$<CONFIG:Debug>:d>")

    # 动态库
    find_package(fmt CONFIG REQUIRED)

    if(DLL)
        #  gmp
        # find_package(PkgConfig REQUIRED)
        # pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)
        list(APPEND EXTRA_LIB fmt::fmt  "${CMAKE_PREFIX_PATH}/$<$<CONFIG:Debug>:debug/>lib/gmp.lib")
    else()
        list(APPEND EXTRA_LIB fmt::fmt mpir)
    endif()

    list(APPEND EXTRA_SRC "${PROJECT_SOURCE_DIR}/src/rand_win.c")
  endif()

elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  message("unix")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Android")
  list(APPEND EXTRA_LIB_PATH "${MY_DEV_LIB_PATH}/android-lib/${ANDROID_ABI}")
  list(APPEND EXTRA_LIB curl)
endif()
message("EXTRA_COMPILE_FLAGS: ${EXTRA_COMPILE_FLAGS}")
message("EXTRA_SRC: ${EXTRA_SRC}")

add_library(common INTERFACE)

# target_link_options(${PROJECT_NAME} PUBLIC ${EXTRA_LINK_FLAGS})
# 编译选项
target_compile_options(common INTERFACE ${EXTRA_COMPILE_FLAGS})
# 宏定义
target_compile_definitions(common INTERFACE ${EXTRA_DEFINES})
# 工程下的头文件的路径
target_include_directories(common INTERFACE ${EXTRA_INCLUDE})
# 第三方库的路径
target_link_directories(common INTERFACE ${EXTRA_LIB_PATH})
# 第三方库
target_link_libraries(common INTERFACE ${EXTRA_LIB})

add_executable(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} PRIVATE common)
# 源文件
list(APPEND EXTRA_SRC "${PROJECT_SOURCE_DIR}/test/test_shamir.cpp")
target_sources(${PROJECT_NAME} PUBLIC ${EXTRA_SRC})

# copy need multi dll to exe dir 
add_custom_command(TARGET shamir POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_PREFIX_PATH}/$<$<CONFIG:Debug>:debug/>bin/fmt$<$<CONFIG:Debug>:d>.dll"
    "${CMAKE_PREFIX_PATH}/$<$<CONFIG:Debug>:debug/>bin/gmp-10.dll"
    $<TARGET_FILE_DIR:shamir>)

if(1)
  message("test")
endif()

include(CTest)