



class Node {

    static nodeID = 0;

    constructor(value) {
        this.value = value;
        this.id = Node.nodeID;
        Node.nodeID++;
        this.next = null;   // horizontal child
        this.child = null;  // vertical child
    }

    setNext(next) {
        this.next = next;
    }

    setChild(child) {
        this.child = child;
    }

    hasNext() {
        return this.next !== null;
    }

    toString() {
        return `Node(${this.value})`;
    }

}


class LinkedList {

    constructor(valuesArray=null) {
        this.head = null;
        this.size = 0;

        if (valuesArray !== null) {
            for (let value of valuesArray) {
                this.push(value);
            }
        }
    }

    push(value) {
        const newNode = new Node(value);
        
        if (this.head === null) {
            this.head = newNode;
        } else {
            let node = this.head;
            while (node.hasNext()) {
                node = node.next;
            }
            node.setNext(newNode);
        }
        this.size++;
        return newNode;
    }

    remove(value) {
        /*
        Remove value from the linked list.
        Return boolean indicating if successfully removed
        */
        if (this.head === null) {
            return false;
        } else {
            if (this.head.value === value) {
                this.head = this.head.next;
                this.size--;
                return true;
            } else {
                let node = this.head;
                while (node.hasNext()) {
                    if (node.next.value === value) {
                        node.next = node.next.next;
                        this.size--;
                        return true;
                    }
                    node = node.next;
                }
                return false;
            }
        }
    }

    search(value) {
        let node = this.head;
        while (node !== null) {
            if (node.value === value) {
                return node;
            }
            node = node.next;
        }
        return null;
    }

    toString() {        
        if (this.head === null) {
            return "LinkedList([])";
        } else {
            let result = "[";
            let node = this.head;
            while (node !== null) {
                result += node.toString() + ", ";
                node = node.next;
            }
            result = result.slice(0, -2) + "]";
            return result;
        }
    }

}


class Trie {

    static TERMINATOR = "<end>";

    constructor(keysArray=null) {
        this.root = null;
        this.maxDepth = 0;
        this.size = 0;

        if (keysArray !== null) {
            for (let key of keysArray) {
                this.insertKey(key);
            }
        }
    }

    insertKey(key) {
        if (this.root === null) {
            this.root = new LinkedList();
        }

        this.maxDepth = Math.max(this.maxDepth, key.length);

        let list = this.root;
        for (let i=0; i<key.length; i++) {
            const character = key[i];
            const node = list.search(character);
            if (node === null) {
                const newNode = list.push(character);
                newNode.setChild(new LinkedList());
                list = newNode.child;
            } else {
                list = node.child;
            }

            if (i === key.length - 1 && list.search(Trie.TERMINATOR) === null) {
                list.push(Trie.TERMINATOR);
                this.size++;
            }
        }
    }

    containsKey(key) {
        if (this.root === null) return false;

        let list = this.root;
        for (let character of key) {
            const node = list.search(character);
            if (node === null) return false;
            list = node.child;
        }
        return list.search(Trie.TERMINATOR) !== null;
    }

    containsPrefix(prefix) {
        if (this.root === null) return false;

        let list = this.root;
        for (let character of prefix) {
            const node = list.search(character);
            if (node === null) return false;
            list = node.child;
        }
        return true;
    }

    remove(key) {
        /*
        Remove key from trie, return boolean indicating if the removal was successful
        Note: does NOT remove prefixes.
        For example:
        const trie = new Trie(["dog"]);
        trie.remove("dog"); // true
        trie.containsKey("dog"); // false
        trie.containsPrefix("dog"); // true
        */
        if (this.root === null) return false;
        let list = this.root;
        for (let character of key) {
            const node = list.search(character);
            if (node === null) return false;
            list = node.child;
        }
        const removalResult = list.remove(Trie.TERMINATOR);
        if (removalResult) this.size--;
        return removalResult;
    }

    prefixSearch(prefix, maxResults=100) {
        /**
         * Search the trie for keys with the given prefix.
         * Returns no more than maxResults keys
         */
        if (this.root === null) return [];

        let list = this.root, node;
        for (let character of prefix) {
            node = list.search(character);
            if (node === null) return []; // prefix not in trie
            list = node.child;
        }

        const results = [];
        if (list.search(Trie.TERMINATOR) !== null) results.push(prefix); 

        // depth-first search for matching keys
        const queue = [{node: list.head, suffix: ""}];
        let current;
        while (results.length < maxResults && queue.length > 0) {
            current = queue.pop();
            let currentNode = current.node, currentSuffix = current.suffix;
            while (currentNode !== null) {
                if (currentNode.value === Trie.TERMINATOR) {
                    currentNode = currentNode.next;
                    continue;
                }
                const searchResult = currentNode.child.search(Trie.TERMINATOR);
                if (searchResult !== null) results.push(prefix + currentSuffix + currentNode.value);
                queue.unshift({
                    node: currentNode.child.head,
                    suffix: currentSuffix + currentNode.value,
                });
                currentNode = currentNode.next;
                if (results.length === maxResults) break;
            }
        }
        return results;
    }

}