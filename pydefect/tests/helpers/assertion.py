# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import json

from monty.json import MSONable, MontyDecoder


def assert_msonable(obj):
    assert isinstance(obj, MSONable)
    assert obj.as_dict() == obj.__class__.from_dict(obj.as_dict()).as_dict()
    assert json.loads(obj.to_json(), cls=MontyDecoder)