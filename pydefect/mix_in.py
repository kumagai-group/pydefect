# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path


class ToJsonFileMixIn:
    def to_json_file(self, filename: str) -> None:
#        def to_json_file(self, filename: Optional[str] = None) -> None:
        Path(filename).write_text(self.to_json())

    # def _filename(self):
    #     self.__class__