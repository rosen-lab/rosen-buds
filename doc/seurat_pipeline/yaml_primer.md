## A brief introduction to YAML
The pipeline's only input is a parameter file encoded in YAML format. YAML is a standardized data format, and many tutorials for it can be found online, such as [this one](http://docs.ansible.com/ansible/latest/YAMLSyntax.html). Processing of YAML files is performed using the R package [yaml](https://github.com/viking/r-yaml). Most text editors with built in syntax highlighting will support YAML files. YAML consists of two basic data structures, lists and key-value maps, which are explained briefly below.

### Lists
Lists simply store multiple pieces of data together. It consists of a name for the list followed by a colon, and then on individual lines the elements of the list preceded by a dash. Additionally, the elements must be indented using spaces (not tabs!). For example, a list named groceries containing elements milk, eggs, and cheese would be written as:
```yaml
groceries:
  - milk
  - eggs
  - cheese
```

There is also a shorthand way to write a list on a single line. This syntax starts the same with the list name followed by a colon. After the colon, the list elements are put in square brackets and separated by commas. So we could rewrite our grocery list as:
```yaml
groceries: [milk, eggs, cheese]
```

### Maps
Maps are similar to lists, but allow you to store a key name as well as a value. The header is similar with a name for the map followed by a colon. On new lines, simply put the key name followed by a colon and then the value. Once again, new lines must be indented using spaces. So if we wanted to encode 1 gallon of milk, 1 dozen eggs, and 1 block of cheese in our grocery list, we would write it as:
```yaml
groceries:
  milk: 1 gallon
  eggs: 1 dozen
  cheese: 1 block
```

As with lists, there is also a shorthand way to write maps. The key-value pairs are put with in curly braces and comma separated. In this syntax, our grocery list looks like:
```yaml
groceries: {milk: 1 gallon, eggs: 1 dozen, cheese: 1 block}
```

### Nesting
Lists and maps can be nested within themselves and each other. So you can put a list in a list, a map in a map, a list in a map, and a map in a list. In truth, every list is actually just a list within a map with the list name as the key. Here are quick examples of the other 3 cases using both long and short form syntax.

### Lists in lists
Long form:
```yaml
list:
  -
    - item 1
    - item 2
  -
    - item 3
    - item 4
```
Short form:
```yaml
list: [[item 1, item 2], [item 3, item 4]]
```

### Maps in maps
Long form:
```yaml
map:
  map1: 
    key1: item1
    key2: item2
  map2:
    key3: item3
    key4: item4
```
Short form:
```yaml
map: {map1: {key1: item1, key2: item2}, map2: {key3: item3, key4: item4}}
```

### Maps in lists
Long form:
```yaml
list:
  - 
    key1: item1
    key2: item2
  - 
    key3: item3
    key4: item4
```
Short form:
```yaml
list: [{key1: item1, key2: item2}, {key3: item3, key4: item4}]
```

### Special characters
The following are special characters in YAML:
```
: { } [ ] , & * # ? | - < > = ! % @ ` ' " \
```

The # character is used for comments. Any text following a # on a line is ignored.

If you wish to include special characters within a string, you must enclose the string with either single quotes or enclose it with double quotes and escape the special character with `\`. Using double quotes also allows you to enter characters like `\t` or `\n`. To use a quotation mark within the same type of quotation marks, you must write the quotation mark twice, for example:
```
"She said ""Hi!"""
```
