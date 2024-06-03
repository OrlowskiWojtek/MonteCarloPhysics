#!/bin/bash

for file in *; do
  # Sprawdź, czy to plik
  if [ -f "$file" ]; then
    # Usuń ostatnią linię w pliku
    sed -i '$ d' "$file"
  fi
done
