
#[inline]
pub fn vbyte_encode(mut val: u32, out: &mut Vec<u8>) {
    loop {
        if val < 0x80 {
            out.push(val as u8);
            return;
        }
        out.push((val as u8 & 0x7F) | 0x80);
        val >>= 7;
    }
}

#[inline]
pub fn vbyte_decode(data: &[u8], pos: &mut usize) -> u32 {
    let mut val: u32 = 0;
    let mut shift = 0;
    loop {
        let b = data[*pos];
        *pos += 1;
        val |= ((b & 0x7F) as u32) << shift;
        if b & 0x80 == 0 {
            return val;
        }
        shift += 7;
    }
}

pub fn delta_vbyte_encode(sorted_ids: &[u32]) -> Vec<u8> {
    let mut out = Vec::new();
    vbyte_encode(sorted_ids.len() as u32, &mut out);
    let mut prev = 0u32;
    for &id in sorted_ids {
        vbyte_encode(id - prev, &mut out);
        prev = id;
    }
    out
}

pub fn delta_vbyte_decode(data: &[u8], pos: &mut usize) -> Vec<u32> {
    let count = vbyte_decode(data, pos) as usize;
    let mut result = Vec::with_capacity(count);
    let mut prev = 0u32;
    for _ in 0..count {
        let delta = vbyte_decode(data, pos);
        prev += delta;
        result.push(prev);
    }
    result
}

#[inline]
pub fn delta_vbyte_skip(data: &[u8], pos: &mut usize) {
    let count = vbyte_decode(data, pos) as usize;
    for _ in 0..count {
        while data[*pos] & 0x80 != 0 {
            *pos += 1;
        }
        *pos += 1;
    }
}