# Definicje zmiennych
CC = g++
CFLAGS = -Wall
INC_DIR = ./include
SRC_DIR = ./src
OBJ_DIR = ./obj

# Pliki źródłowe i nagłówkowe
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
INCS = $(wildcard $(INC_DIR)/*.h)

# Nazwa pliku wyjściowego
EXEC = exec

# Reguła kompilacji i linkowania
$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@

# Reguła kompilacji plików źródłowych do plików obiektowych
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCS)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

# Czyszczenie plików tymczasowych i wyjściowych
clean:
	rm -rf $(OBJ_DIR)/*.o $(EXEC)

# Reguła "phony" do uruchamiania programu
run: $(EXEC)
	./$(EXEC)

# Domyślna reguła
.PHONY: all
all: $(EXEC)

# Wymagania dla reguł
$(OBJS): | $(OBJ_DIR)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Ustawienie domyślnej reguły
.DEFAULT_GOAL := all