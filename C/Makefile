TARGET_EXEC_NAME ?= exec

CC = gcc
CXX = g++

BUILD_DIR ?= .build
SRC_DIRS ?= src


ifeq ($(OS),Windows_NT)
TARGET_EXEC_NAME ?= $(TARGET_EXEC_NAME).exe
endif

TARGET_EXEC ?= $(BUILD_DIR)/$(TARGET_EXEC_NAME)

ifeq ($(OS),Windows_NT)
RUN ?= $(TARGET)
else
RUN = ./$(TARGET)
endif

SRCS := $(shell find $(SRC_DIRS) -maxdepth 1 -name "*.cpp" -or -name "*.c" ! -name "test.c" -or -name "*.s")
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)


INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP -g -fPIE -Wno-unused-label

LDFLAGS := -lm

$(TARGET_EXEC): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

# assembly
$(BUILD_DIR)/%.s.o: %.s
	$(MKDIR_P) $(dir $@)
	$(AS) $(ASFLAGS) -c $< -o $@

# c source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


.PHONY: clean
.PHONY: run

clean:
	$(RM) -r $(BUILD_DIR)

run:
	./$(TARGET_EXEC)


-include $(DEPS)

MKDIR_P ?= mkdir -p
