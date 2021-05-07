#!/bin/bash

clang-format -i -style=file            \
  src_lib/*                            \
  include_lib/*                        \
  script/*                             \
  BRparity/script/example_brparity.cpp
