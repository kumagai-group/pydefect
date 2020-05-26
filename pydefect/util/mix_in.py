# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import re
from pathlib import Path
from typing import Optional


class ToJsonFileMixIn:
    def to_json_file(self, filename: Optional[str] = None) -> None:
        filename = filename or self._filename
        print(filename)
        Path(filename).write_text(self.to_json())

    @property
    def _filename(self):
        """ ClassForThis -> class_for_this.json"""
        class_name = self.__class__.__name__
        return re.sub("(?<!^)(?=[A-Z])", "_", class_name).lower() + ".json"
