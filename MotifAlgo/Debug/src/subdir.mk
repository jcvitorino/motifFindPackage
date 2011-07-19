################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/OligoAnalisys.cpp \
../src/XNA.cpp \
../src/XNAseq.cpp \
../src/function.cpp \
../src/main.cpp \
../src/medianString.cpp 

OBJS += \
./src/OligoAnalisys.o \
./src/XNA.o \
./src/XNAseq.o \
./src/function.o \
./src/main.o \
./src/medianString.o 

CPP_DEPS += \
./src/OligoAnalisys.d \
./src/XNA.d \
./src/XNAseq.d \
./src/function.d \
./src/main.d \
./src/medianString.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


